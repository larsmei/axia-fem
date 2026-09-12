use std::collections::HashMap;

use crate::beam;
use crate::constraint::{self, DofMap};
use crate::dat;
use crate::elem::{
    element_ke, element_nodal_stress, hex8_body_force, hex8_face_pressure, quad4_body_force,
    quad4_edge_pressure, tet4_body_force, von_mises,
};
use crate::error::{err, Result};
use crate::extra;
use crate::frd;
use crate::heat;
use crate::linalg::solve_kff;
use crate::model::{Dload, ElemKind, FluxKind, InitKind, Model, Procedure};
use crate::quadratic;
use crate::shell;

pub struct SolveOutput {
    pub model: Model,
    pub u: Vec<[f64; 3]>,
    pub ur: Vec<[f64; 3]>,
    pub rf: Vec<[f64; 3]>,
    pub rm: Vec<[f64; 3]>,
    pub stress: Vec<[f64; 6]>,
    pub strain: Vec<[f64; 6]>,
    pub von_mises: Vec<f64>,
    pub stress_gp: Vec<(i32, usize, [f64; 6])>,
    pub frd: String,
    pub dat: String,
    pub nfree: usize,
    pub ndof: usize,
    pub solver: String,
    pub iters: usize,
    pub residual: f64,
    pub time_ms: f64,
    pub procedure: String,
    pub frequencies: Vec<f64>,
    pub buckles: Vec<f64>,
}

fn elem_xyz(model: &Model, nodes: &[i32]) -> Result<Vec<[f64; 3]>> {
    let mut xyz = Vec::with_capacity(nodes.len());
    for &id in nodes {
        let i = model.node_index(id)?;
        xyz.push(model.coords[i]);
    }
    Ok(xyz)
}

fn dof_of(ndn: usize, node_index: usize, dir: usize) -> usize {
    ndn * node_index + dir
}

pub fn solve(model: Model) -> Result<SolveOutput> {
    let t0 = now_ms();
    if matches!(model.procedure, Procedure::HeatTransfer { .. }) {
        return solve_heat(model, t0);
    }
    let nlgeom = matches!(
        model.procedure,
        Procedure::Static { nlgeom: true, .. }
    );
    let truss2_only = !model.elements.is_empty()
        && model.elements.iter().all(|e| e.kind == ElemKind::Truss2);
    if nlgeom || model.has_plastic() {
        if !truss2_only {
            return err(
                "NLGEOM und *PLASTIC sind in 1.0 nur für T3D2-Fachwerke implementiert.",
            );
        }
        return solve_truss_newton(model, t0, nlgeom);
    }
    solve_linear(model, t0)
}

fn solve_linear(model: Model, t0: f64) -> Result<SolveOutput> {
    let ndn = model.ndof_node();
    let nnode = model.node_ids.len();
    let ndof = ndn * nnode;
    if ndof == 0 {
        return err("Modell ohne Freiheitsgrade.");
    }

    let mut struct_node = vec![false; nnode];
    if ndn == 6 {
        for el in &model.elements {
            if el.kind.is_beam() || el.kind.is_shell() {
                for &id in &el.nodes {
                    struct_node[model.node_index(id)?] = true;
                }
            }
        }
        for rb in &model.rigid_bodies {
            if let Ok(i) = model.node_index(rb.ref_node) {
                struct_node[i] = true;
            }
        }
        for c in &model.couplings {
            if c.kinematic {
                if let Ok(i) = model.node_index(c.ref_node) {
                    struct_node[i] = true;
                }
            }
        }
    }

    // Constrained dofs: last BC wins
    let mut prescribed: HashMap<usize, f64> = HashMap::new();
    for bc in &model.bcs {
        if bc.dof >= ndn {
            continue;
        }
        let ni = model.node_index(bc.node)?;
        prescribed.insert(dof_of(ndn, ni, bc.dof), bc.value);
    }
    // Unused rotational DOFs on continuum-only nodes would be singular.
    if ndn == 6 {
        for ni in 0..nnode {
            if !struct_node[ni] {
                for r in 3..6 {
                    prescribed.entry(dof_of(ndn, ni, r)).or_insert(0.0);
                }
            }
        }
    }

    let mut f_full = vec![0.0; ndof];
    let mut m_full = vec![0.0; ndof];

    let mut trips: Vec<(usize, usize, f64)> = Vec::new();
    for el in &model.elements {
        let xyz = elem_xyz(&model, &el.nodes)?;
        let mat = model.material_for(el)?;
        let th = if el.kind.is_spring() {
            model.spring_k_for(el)?
        } else {
            model.thickness_for(el)
        };
        let sec = if el.kind.is_beam() {
            Some(model.beam_section_for(el)?)
        } else {
            None
        };
        let mut kef = element_ke(el.kind, &xyz, mat.e, mat.nu, th, sec.as_ref())?;
        if !model.node_transform.is_empty() {
            constraint::transform_ke(
                &mut kef.ke,
                el.kind.nnodes(),
                el.kind.ndof_per_node(),
                &el.nodes,
                &model.node_transform,
            );
        }
        let nn = el.kind.nnodes();
        let local_dim = el.kind.ndof_per_node();
        let mut gdofs = Vec::with_capacity(nn * local_dim);
        for a in 0..nn {
            let ni = model.node_index(el.nodes[a])?;
            for d in 0..local_dim {
                gdofs.push(dof_of(ndn, ni, d));
            }
        }
        let m = gdofs.len();
        for i in 0..m {
            for j in 0..m {
                let v = kef.ke[i * kef.ndof + j];
                if v.abs() > 0.0 {
                    trips.push((gdofs[i], gdofs[j], v));
                }
            }
        }
        if mat.density.abs() > 0.0 && kef.volume.abs() > 0.0 {
            let mnode = mat.density * kef.volume / nn as f64;
            for a in 0..nn {
                let ni = model.node_index(el.nodes[a])?;
                for d in 0..3.min(local_dim) {
                    m_full[dof_of(ndn, ni, d)] += mnode;
                }
                if local_dim >= 6 {
                    let c2 = kef.volume.abs().powf(2.0 / 3.0).max(1e-6);
                    for r in 3..6 {
                        m_full[dof_of(ndn, ni, r)] += mnode * c2 * 1e-6;
                    }
                }
            }
        }

        for dl in &model.dloads {
            match dl {
                Dload::Pressure { elem, face, mag } if *elem == el.id => {
                    if el.kind.is_beam() {
                        if let Some(sec) = sec.as_ref() {
                            // P1/P2/P3 on a beam → local distributed load
                            if *face >= 1 {
                                let fe = beam::line_load_local(
                                    el.kind,
                                    &xyz,
                                    sec,
                                    (*face as usize).saturating_sub(1),
                                    *mag,
                                )?;
                                scatter_fe(&fe, &gdofs, 6, local_dim, &mut f_full);
                            }
                        }
                    } else if el.kind.is_shell() {
                        let fe = shell::pressure_force(el.kind, &xyz, *mag)?;
                        scatter_fe(&fe, &gdofs, 6, local_dim, &mut f_full);
                    } else {
                        apply_pressure(
                            &model,
                            el.kind,
                            &xyz,
                            *face,
                            *mag,
                            th,
                            &gdofs,
                            local_dim,
                            &mut f_full,
                        )?;
                    }
                }
                Dload::Grav { mag, dir } => {
                    let mut ndir = *dir;
                    let len = (ndir[0] * ndir[0] + ndir[1] * ndir[1] + ndir[2] * ndir[2]).sqrt();
                    if len > 0.0 {
                        ndir[0] /= len;
                        ndir[1] /= len;
                        ndir[2] /= len;
                    }
                    let bx = mat.density * *mag * ndir[0];
                    let by = mat.density * *mag * ndir[1];
                    let bz = mat.density * *mag * ndir[2];
                    if el.kind.is_beam() {
                        if let Some(sec) = sec.as_ref() {
                            let fe = beam::body_force(el.kind, &xyz, sec, bx, by, bz)?;
                            scatter_fe(&fe, &gdofs, 6, local_dim, &mut f_full);
                        }
                    } else if el.kind.is_shell() {
                        let fe = shell::body_force(el.kind, &xyz, bx, by, bz, th)?;
                        scatter_fe(&fe, &gdofs, 6, local_dim, &mut f_full);
                    } else {
                        apply_body(
                            &model,
                            el.kind,
                            &xyz,
                            bx,
                            by,
                            bz,
                            th,
                            &gdofs,
                            local_dim,
                            &mut f_full,
                        )?;
                    }
                }
                Dload::BeamGlobal { elem, dir, mag } if *elem == el.id && el.kind.is_beam() => {
                    let fe = beam::line_load_global(
                        el.kind,
                        &xyz,
                        mag * dir[0],
                        mag * dir[1],
                        mag * dir[2],
                    )?;
                    scatter_fe(&fe, &gdofs, 6, local_dim, &mut f_full);
                }
                _ => {}
            }
        }
        if mat.alpha.abs() > 0.0 {
            let mut tsum = 0.0;
            for &id in &el.nodes {
                tsum += model.temperature_at(id);
            }
            let dt_th = tsum / nn as f64 - mat.tref;
            if dt_th.abs() > 0.0 {
                if el.kind.is_truss() {
                    let fe = extra::truss_thermal_force(&xyz, mat.e, th, mat.alpha, dt_th);
                    scatter_fe(&fe, &gdofs, 3, local_dim, &mut f_full);
                } else if matches!(
                    el.kind,
                    ElemKind::Hex8 | ElemKind::Hex8I | ElemKind::Hex8R
                ) {
                    let fe = extra::hex8_thermal_force(&xyz, mat.e, mat.nu, mat.alpha, dt_th)?;
                    scatter_fe(&fe, &gdofs, 3, local_dim, &mut f_full);
                }
            }
        }
    }

    add_cloads(&model, ndn, 0.0, &mut f_full)?;
    let f_dload = {
        let mut f = f_full.clone();
        // strip t=0 cloads so Dynamic can re-apply with amplitude(t)
        let mut clo = vec![0.0; ndof];
        add_cloads(&model, ndn, 0.0, &mut clo)?;
        for i in 0..ndof {
            f[i] -= clo[i];
        }
        f
    };

    let mpcs = constraint::build_all_mpcs(&model, ndn)?;
    let map = DofMap::build(ndof, &prescribed, &mpcs)?;
    let nfree = map.n_ind;
    let mut solver = "prescribed".to_string();
    let mut iters = 0usize;
    let mut residual = 0.0;
    let mut frequencies = Vec::new();
    let mut buckles = Vec::new();

    let mut u_full = if nfree == 0 {
        map.u0.clone()
    } else if let Procedure::Frequency { nmodes } = model.procedure {
        let mut m_ind = vec![0.0; nfree];
        for d in 0..ndof {
            for &(a, ta) in &map.t_row[d] {
                m_ind[a] += ta * ta * m_full[d];
            }
        }
        if m_ind.iter().all(|v| *v <= 0.0) {
            return err("*FREQUENCY: *DENSITY fehlt oder Masse ist null.");
        }
        let (ff_trips, _) = map.reduce(&trips, &f_full);
        let ev = crate::eigen::subspace_gen(nfree, ff_trips, &m_ind, nmodes.max(1))?;
        frequencies = ev
            .values
            .iter()
            .map(|l| {
                if *l > 0.0 {
                    l.sqrt() / (2.0 * std::f64::consts::PI)
                } else {
                    0.0
                }
            })
            .collect();
        solver = format!("eigen ({} modes, {})", frequencies.len(), "subspace");
        iters = 1;
        residual = 0.0;
        let mode0 = ev.vectors.first().cloned().unwrap_or_else(|| vec![0.0; nfree]);
        map.reconstruct(&mode0)
    } else if let Procedure::Dynamic { dt, period } = model.procedure {
        let (u_dyn, name, it, res) = newmark(
            &model,
            ndn,
            ndof,
            nfree,
            &map,
            &trips,
            &m_full,
            &f_dload,
            dt,
            period,
        )?;
        solver = name;
        iters = it;
        residual = res;
        u_dyn
    } else {
        let (ff_trips, rhs) = map.reduce(&trips, &f_full);
        let solved = solve_kff(nfree, ff_trips, &rhs)?;
        solver = solved.name;
        iters = solved.iters;
        residual = solved.residual;
        map.reconstruct(&solved.x)
    };

    // Reactions in the (possibly local) analysis DOFs, then rotate to global.
    let mut ku = vec![0.0; ndof];
    for (i, j, v) in &trips {
        ku[*i] += *v * u_full[*j];
    }
    let mut rf_full = vec![0.0; ndof];
    for d in 0..ndof {
        rf_full[d] = ku[d] - f_full[d];
    }
    constraint::dofs_to_global(&mut u_full, ndn, &model.node_ids, &model.node_transform);
    constraint::dofs_to_global(&mut rf_full, ndn, &model.node_ids, &model.node_transform);

    let mut u = vec![[0.0; 3]; nnode];
    let mut ur = vec![[0.0; 3]; nnode];
    let mut rf = vec![[0.0; 3]; nnode];
    let mut rm = vec![[0.0; 3]; nnode];
    for ni in 0..nnode {
        for d in 0..3.min(ndn) {
            u[ni][d] = u_full[dof_of(ndn, ni, d)];
            rf[ni][d] = rf_full[dof_of(ndn, ni, d)];
        }
        if ndn >= 6 {
            for d in 0..3 {
                ur[ni][d] = u_full[dof_of(ndn, ni, 3 + d)];
                rm[ni][d] = rf_full[dof_of(ndn, ni, 3 + d)];
            }
        }
    }

    // Nodal averaged stresses
    let mut acc = vec![[0.0; 6]; nnode];
    let mut cnt = vec![0.0; nnode];
    let mut stress_gp = Vec::new();
    for el in &model.elements {
        let xyz = elem_xyz(&model, &el.nodes)?;
        let mat = model.material_for(el)?;
        let nn = el.kind.nnodes();
        let local_dim = el.kind.ndof_per_node();
        let th = model.thickness_for(el);
        let sec = if el.kind.is_beam() {
            Some(model.beam_section_for(el)?)
        } else {
            None
        };
        let mut ue = vec![0.0; nn * local_dim];
        for a in 0..nn {
            let ni = model.node_index(el.nodes[a])?;
            for d in 0..local_dim {
                ue[a * local_dim + d] = u_full[dof_of(ndn, ni, d)];
            }
        }
        let sn = element_nodal_stress(el.kind, &xyz, &ue, mat.e, mat.nu, sec.as_ref(), th)?;
        for a in 0..nn {
            let ni = model.node_index(el.nodes[a])?;
            for c in 0..6 {
                acc[ni][c] += sn[a][c];
            }
            cnt[ni] += 1.0;
        }
        let mut mean = [0.0; 6];
        for a in 0..nn {
            for c in 0..6 {
                mean[c] += sn[a][c] / nn as f64;
            }
        }
        stress_gp.push((el.id, 1usize, mean));
    }
    let mut stress = vec![[0.0; 6]; nnode];
    let mut strain = vec![[0.0; 6]; nnode];
    let mut vm = vec![0.0; nnode];
    for i in 0..nnode {
        if cnt[i] > 0.0 {
            for c in 0..6 {
                stress[i][c] = acc[i][c] / cnt[i];
            }
        }
        vm[i] = von_mises(&stress[i]);
        let mat = model.materials.values().next().copied().unwrap_or_default();
        let e = mat.e;
        let nu = mat.nu;
        let tr = stress[i][0] + stress[i][1] + stress[i][2];
        strain[i][0] = ((1.0 + nu) * stress[i][0] - nu * tr) / e;
        strain[i][1] = ((1.0 + nu) * stress[i][1] - nu * tr) / e;
        strain[i][2] = ((1.0 + nu) * stress[i][2] - nu * tr) / e;
        let g2 = e / (1.0 + nu);
        strain[i][3] = stress[i][3] / g2;
        strain[i][4] = stress[i][4] / g2;
        strain[i][5] = stress[i][5] / g2;
    }

    if let Procedure::Buckle { nmodes } = model.procedure {
        if nfree > 0 {
            let mut kg_trips: Vec<(usize, usize, f64)> = Vec::new();
            for el in &model.elements {
                let xyz = elem_xyz(&model, &el.nodes)?;
                let nn = el.kind.nnodes();
                let local_dim = el.kind.ndof_per_node();
                let mut ue = vec![0.0; nn * local_dim];
                let mut gdofs = Vec::new();
                for a in 0..nn {
                    let ni = model.node_index(el.nodes[a])?;
                    for d in 0..local_dim {
                        gdofs.push(dof_of(ndn, ni, d));
                        ue[a * local_dim + d] = u_full[dof_of(ndn, ni, d)];
                    }
                }
                let kg = if el.kind.is_truss() {
                    let area = model.thickness_for(el);
                    let mat = model.material_for(el)?;
                    let i1 = if nn == 2 { 1 } else { 2 };
                    let mut d = [
                        xyz[i1][0] - xyz[0][0],
                        xyz[i1][1] - xyz[0][1],
                        xyz[i1][2] - xyz[0][2],
                    ];
                    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt().max(1e-18);
                    d[0] /= len;
                    d[1] /= len;
                    d[2] /= len;
                    let du = [
                        ue[3 * i1] - ue[0],
                        ue[3 * i1 + 1] - ue[1],
                        ue[3 * i1 + 2] - ue[2],
                    ];
                    let n_ax = mat.e * area * (du[0] * d[0] + du[1] * d[1] + du[2] * d[2]) / len;
                    crate::eigen::truss_kg(&xyz, n_ax)
                } else if el.kind.is_beam() {
                    let sec = model.beam_section_for(el)?;
                    let mat = model.material_for(el)?;
                    let dx = xyz[1][0] - xyz[0][0];
                    let dy = xyz[1][1] - xyz[0][1];
                    let dz = xyz[1][2] - xyz[0][2];
                    let len = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-18);
                    let axial = (ue[6] - ue[0]) * dx / len
                        + (ue[7] - ue[1]) * dy / len
                        + (ue[8] - ue[2]) * dz / len;
                    let n_ax = mat.e * sec.area * axial / len;
                    crate::eigen::beam_kg(&xyz, n_ax, nn)?
                } else if matches!(el.kind, ElemKind::Hex8 | ElemKind::Hex8I | ElemKind::Hex8R) {
                    let mut mean = [0.0; 6];
                    if let Some((_, _, s)) = stress_gp.iter().find(|(id, _, _)| *id == el.id) {
                        mean = *s;
                    }
                    crate::eigen::hex8_kg(&xyz, &mean)?
                } else {
                    continue;
                };
                let m = gdofs.len();
                for i in 0..m {
                    for j in 0..m {
                        let v = kg[i * m + j];
                        if v.abs() > 0.0 {
                            kg_trips.push((gdofs[i], gdofs[j], v));
                        }
                    }
                }
            }
            let (k_ff, _) = map.reduce(&trips, &f_full);
            let (kg_ff, _) = map.reduce(&kg_trips, &f_full);
            let a_ff: Vec<_> = kg_ff.into_iter().map(|(i, j, v)| (i, j, -v)).collect();
            match crate::eigen::subspace_ab(nfree, k_ff, a_ff, nmodes.max(1)) {
                Ok(ev) => {
                    buckles = ev.values;
                    solver = format!("buckle ({} factors)", buckles.len());
                    if let Some(v0) = ev.vectors.first() {
                        u_full = map.reconstruct(v0);
                    }
                }
                Err(e) => return Err(e),
            }
            for ni in 0..nnode {
                for d in 0..3.min(ndn) {
                    u[ni][d] = u_full[dof_of(ndn, ni, d)];
                }
                if ndn >= 6 {
                    for d in 0..3 {
                        ur[ni][d] = u_full[dof_of(ndn, ni, 3 + d)];
                    }
                }
            }
        }
    }

    let frd_s = frd::write_frd(&model, &u, &stress, &rf, &strain);
    let dat_s = dat::write_dat(&model, &u, &stress_gp, &rf);
    let dt = now_ms() - t0;
    let procedure = model.procedure.name().to_string();

    Ok(SolveOutput {
        model,
        u,
        ur,
        rf,
        rm,
        stress,
        strain,
        von_mises: vm,
        stress_gp,
        frd: frd_s,
        dat: dat_s,
        nfree,
        ndof,
        solver,
        iters,
        residual,
        time_ms: dt,
        procedure,
        frequencies,
        buckles,
    })
}

fn apply_pressure(
    _model: &Model,
    kind: ElemKind,
    xyz: &[[f64; 3]],
    face: i32,
    mag: f64,
    th: f64,
    gdofs: &[usize],
    local_dim: usize,
    f_full: &mut [f64],
) -> Result<()> {
    match kind {
        ElemKind::Hex8 | ElemKind::Hex8I | ElemKind::Hex8R => {
            let mut p = [[0.0; 3]; 8];
            for i in 0..8 {
                p[i] = xyz[i];
            }
            let fe = hex8_face_pressure(&p, face, mag)?;
            scatter_fe(&fe, gdofs, 3, local_dim, f_full);
        }
        ElemKind::Wedge6 => {
            let fe = extra::wedge6_face_pressure(xyz, face, mag)?;
            scatter_fe(&fe, gdofs, 3, local_dim, f_full);
        }
        ElemKind::Hex20 | ElemKind::Hex20R => {
            let fe = quadratic::hex20_face_pressure(xyz, face, mag)?;
            scatter_fe(&fe, gdofs, 3, local_dim, f_full);
        }
        ElemKind::Quad4Ps | ElemKind::Quad4Pe => {
            if face <= 0 {
                return Ok(());
            }
            let mut p = [[0.0; 2]; 4];
            for i in 0..4 {
                p[i] = [xyz[i][0], xyz[i][1]];
            }
            let fe = quad4_edge_pressure(&p, face, mag, th)?;
            scatter_fe(&fe, gdofs, 2, local_dim, f_full);
        }
        ElemKind::Quad8Ps | ElemKind::Quad8Pe | ElemKind::Quad8RPs | ElemKind::Quad8RPe => {
            if face <= 0 {
                return Ok(());
            }
            let mut p = [[0.0; 2]; 8];
            for i in 0..8 {
                p[i] = [xyz[i][0], xyz[i][1]];
            }
            let fe = quadratic::quad8_edge_pressure(&p, face, mag, th)?;
            scatter_fe(&fe, gdofs, 2, local_dim, f_full);
        }
        _ => {}
    }
    Ok(())
}

fn apply_body(
    _model: &Model,
    kind: ElemKind,
    xyz: &[[f64; 3]],
    bx: f64,
    by: f64,
    bz: f64,
    th: f64,
    gdofs: &[usize],
    local_dim: usize,
    f_full: &mut [f64],
) -> Result<()> {
    match kind {
        ElemKind::Hex8 | ElemKind::Hex8I | ElemKind::Hex8R => {
            let mut p = [[0.0; 3]; 8];
            for i in 0..8 {
                p[i] = xyz[i];
            }
            let fe = hex8_body_force(&p, bx, by, bz)?;
            scatter_fe(&fe, gdofs, 3, local_dim, f_full);
        }
        ElemKind::Wedge6 => {
            let fe = extra::wedge6_body_force(xyz, bx, by, bz)?;
            scatter_fe(&fe, gdofs, 3, local_dim, f_full);
        }
        ElemKind::Truss2 | ElemKind::Truss3 => {
            let fe = extra::truss_body_force(xyz, th, bx, by, bz);
            scatter_fe(&fe, gdofs, 3, local_dim, f_full);
        }
        ElemKind::Tet4 => {
            let mut p = [[0.0; 3]; 4];
            for i in 0..4 {
                p[i] = xyz[i];
            }
            let fe = tet4_body_force(&p, bx, by, bz)?;
            scatter_fe(&fe, gdofs, 3, local_dim, f_full);
        }
        ElemKind::Quad4Ps | ElemKind::Quad4Pe => {
            let mut p = [[0.0; 2]; 4];
            for i in 0..4 {
                p[i] = [xyz[i][0], xyz[i][1]];
            }
            let fe = quad4_body_force(&p, bx, by, th)?;
            scatter_fe(&fe, gdofs, 2, local_dim, f_full);
        }
        ElemKind::Hex20 | ElemKind::Hex20R => {
            let fe = quadratic::hex20_body_force(xyz, bx, by, bz, kind.reduced_int())?;
            scatter_fe(&fe, gdofs, 3, local_dim, f_full);
        }
        ElemKind::Tet10 => {
            let fe = quadratic::tet10_body_force(xyz, bx, by, bz)?;
            scatter_fe(&fe, gdofs, 3, local_dim, f_full);
        }
        ElemKind::Quad8Ps | ElemKind::Quad8Pe | ElemKind::Quad8RPs | ElemKind::Quad8RPe => {
            let mut p = [[0.0; 2]; 8];
            for i in 0..8 {
                p[i] = [xyz[i][0], xyz[i][1]];
            }
            let fe = quadratic::quad8_body_force(&p, bx, by, th, kind.reduced_int())?;
            scatter_fe(&fe, gdofs, 2, local_dim, f_full);
        }
        ElemKind::Tri6Ps | ElemKind::Tri6Pe => {
            let mut p = [[0.0; 2]; 6];
            for i in 0..6 {
                p[i] = [xyz[i][0], xyz[i][1]];
            }
            let fe = quadratic::tri6_body_force(&p, bx, by, th)?;
            scatter_fe(&fe, gdofs, 2, local_dim, f_full);
        }
        _ => {}
    }
    Ok(())
}

fn scatter_fe(fe: &[f64], gdofs: &[usize], fe_dim: usize, local_dim: usize, f_full: &mut [f64]) {
    let nn = gdofs.len() / local_dim;
    for a in 0..nn {
        for d in 0..local_dim.min(fe_dim) {
            f_full[gdofs[a * local_dim + d]] += fe[a * fe_dim + d];
        }
    }
}

fn add_cloads(model: &Model, ndn: usize, t: f64, f: &mut [f64]) -> Result<()> {
    for c in &model.cloads {
        if c.dof >= ndn {
            continue;
        }
        let ni = model.node_index(c.node)?;
        let s = model.amp_value(&c.amplitude, t);
        f[dof_of(ndn, ni, c.dof)] += c.mag * s;
    }
    Ok(())
}

fn newmark(
    model: &Model,
    ndn: usize,
    ndof: usize,
    nfree: usize,
    map: &DofMap,
    trips: &[(usize, usize, f64)],
    m_full: &[f64],
    f_dload: &[f64],
    dt_in: f64,
    period: f64,
) -> Result<(Vec<f64>, String, usize, f64)> {
    if nfree == 0 {
        return Ok((map.u0.clone(), "prescribed".into(), 0, 0.0));
    }
    let mut m_ind = vec![0.0; nfree];
    for d in 0..ndof {
        for &(a, ta) in &map.t_row[d] {
            m_ind[a] += ta * ta * m_full[d];
        }
    }
    if m_ind.iter().all(|v| *v <= 0.0) {
        return err("*DYNAMIC: *DENSITY fehlt oder Masse ist null.");
    }
    let nsteps = ((period / dt_in).round() as usize).clamp(1, 20_000);
    let dt = period / nsteps as f64;
    let beta = 0.25_f64;
    let gamma = 0.5_f64;
    let a0 = 1.0 / (beta * dt * dt);
    let a1 = gamma / (beta * dt);
    let a2 = 1.0 / (beta * dt);
    let a3 = 1.0 / (2.0 * beta) - 1.0;
    let a4 = gamma / beta - 1.0;
    let a5 = dt * (gamma / (2.0 * beta) - 1.0);
    let alpha_r = model.damp_alpha;
    let beta_r = model.damp_beta;

    let k_scale = 1.0 + a1 * beta_r;
    let m_scale = a0 + a1 * alpha_r;
    let mut keff_trips: Vec<(usize, usize, f64)> = trips
        .iter()
        .map(|&(i, j, v)| (i, j, v * k_scale))
        .collect();
    for d in 0..ndof {
        if m_full[d].abs() > 0.0 {
            keff_trips.push((d, d, m_full[d] * m_scale));
        }
    }

    let mut u_full = map.u0.clone();
    for ic in &model.init {
        if ic.kind != InitKind::Displacement || ic.dof >= ndn {
            continue;
        }
        if let Ok(ni) = model.node_index(ic.node) {
            let d = dof_of(ndn, ni, ic.dof);
            if d < ndof && map.ind_of[d] >= 0 {
                u_full[d] = ic.value;
            }
        }
    }
    let mut v_full = vec![0.0; ndof];
    for ic in &model.init {
        if ic.kind != InitKind::Velocity || ic.dof >= ndn {
            continue;
        }
        if let Ok(ni) = model.node_index(ic.node) {
            let d = dof_of(ndn, ni, ic.dof);
            if d < ndof {
                v_full[d] = ic.value;
            }
        }
    }

    let mut f0 = f_dload.to_vec();
    add_cloads(model, ndn, 0.0, &mut f0)?;
    let mut ku = vec![0.0; ndof];
    let mut kv = vec![0.0; ndof];
    for &(i, j, v) in trips {
        ku[i] += v * u_full[j];
        kv[i] += v * v_full[j];
    }
    let mut a_full = vec![0.0; ndof];
    for i in 0..ndof {
        if map.ind_of[i] < 0 {
            v_full[i] = 0.0;
            continue;
        }
        let rhs = f0[i] - ku[i] - alpha_r * m_full[i] * v_full[i] - beta_r * kv[i];
        if m_full[i].abs() > 1e-30 {
            a_full[i] = rhs / m_full[i];
        }
    }

    let mut solver_name = String::new();
    let mut residual = 0.0;
    for step in 1..=nsteps {
        let t = step as f64 * dt;
        let mut f_t = f_dload.to_vec();
        add_cloads(model, ndn, t, &mut f_t)?;
        let mut reff = f_t;
        for i in 0..ndof {
            reff[i] += m_full[i] * (a0 * u_full[i] + a2 * v_full[i] + a3 * a_full[i]);
        }
        let mut pred = vec![0.0; ndof];
        for i in 0..ndof {
            pred[i] = a1 * u_full[i] + a4 * v_full[i] + a5 * a_full[i];
        }
        for i in 0..ndof {
            reff[i] += alpha_r * m_full[i] * pred[i];
        }
        let mut kpred = vec![0.0; ndof];
        for &(i, j, v) in trips {
            kpred[i] += v * pred[j];
        }
        for i in 0..ndof {
            reff[i] += beta_r * kpred[i];
        }
        let (ff, rhs) = map.reduce(&keff_trips, &reff);
        let solved = solve_kff(nfree, ff, &rhs)?;
        solver_name = solved.name;
        residual = solved.residual;
        let u_new = map.reconstruct(&solved.x);
        let mut a_new = vec![0.0; ndof];
        let mut v_new = vec![0.0; ndof];
        for i in 0..ndof {
            if map.ind_of[i] < 0 {
                continue;
            }
            a_new[i] = a0 * (u_new[i] - u_full[i]) - a2 * v_full[i] - a3 * a_full[i];
            v_new[i] = v_full[i] + dt * ((1.0 - gamma) * a_full[i] + gamma * a_new[i]);
        }
        u_full = u_new;
        v_full = v_new;
        a_full = a_new;
    }
    Ok((
        u_full,
        format!("Newmark ({solver_name}, {nsteps} steps)"),
        nsteps,
        residual,
    ))
}

fn solve_heat(model: Model, t0: f64) -> Result<SolveOutput> {
    let nnode = model.node_ids.len();
    let ndof = nnode;
    let mut prescribed: HashMap<usize, f64> = HashMap::new();
    for bc in &model.thermal_bcs {
        let ni = model.node_index(bc.node)?;
        prescribed.insert(ni, bc.value);
    }
    for ic in &model.init {
        if ic.kind == InitKind::Temperature {
            let ni = model.node_index(ic.node)?;
            prescribed.entry(ni).or_insert(ic.value);
        }
    }
    if prescribed.is_empty() && model.films.is_empty() {
        return err("*HEAT TRANSFER: keine Temperatur-Randbedingung (DOF 11) und kein *FILM.");
    }

    let mut trips: Vec<(usize, usize, f64)> = Vec::new();
    let mut f = vec![0.0; ndof];
    let mut c_diag = vec![0.0; ndof];
    for el in &model.elements {
        let xyz = elem_xyz(&model, &el.nodes)?;
        let mat = model.material_for(el)?;
        let kth = mat.conductivity;
        if kth <= 0.0 {
            return err(format!(
                "*HEAT TRANSFER: *CONDUCTIVITY fehlt für Element {} (ELSET={}).",
                el.id, el.elset
            ));
        }
        let area = if el.kind.is_beam() {
            model.beam_section_for(el)?.area
        } else {
            model.thickness_for(el)
        };
        let (ke, vol) = heat::element_conductivity(el.kind, &xyz, kth, area)?;
        let nn = el.kind.nnodes();
        let mut gdofs = Vec::with_capacity(nn);
        for &id in &el.nodes {
            gdofs.push(model.node_index(id)?);
        }
        for i in 0..nn {
            for j in 0..nn {
                let v = ke[i * nn + j];
                if v.abs() > 0.0 {
                    trips.push((gdofs[i], gdofs[j], v));
                }
            }
        }
        if mat.density.abs() > 0.0 && mat.specific_heat.abs() > 0.0 && vol.abs() > 0.0 {
            let ci = mat.density * mat.specific_heat * vol / nn as f64;
            for &g in &gdofs {
                c_diag[g] += ci;
            }
        }
        for df in &model.dfluxes {
            if df.elem != el.id {
                continue;
            }
            match df.kind {
                FluxKind::Body => {
                    let fe = heat::element_body_heat(el.kind, &xyz, df.mag, area)?;
                    for a in 0..nn.min(fe.len()) {
                        f[gdofs[a]] += fe[a];
                    }
                }
                FluxKind::Face(face) if face >= 1 => {
                    if matches!(
                        el.kind,
                        ElemKind::Hex8 | ElemKind::Hex8I | ElemKind::Hex8R
                    ) {
                        let (_ke_f, fe) = heat::hex8_face_heat(&xyz, face, df.mag, 0.0, 0.0)?;
                        for a in 0..8 {
                            f[gdofs[a]] += fe[a];
                        }
                    }
                }
                _ => {}
            }
        }
        for fm in &model.films {
            if fm.elem != el.id {
                continue;
            }
            if matches!(
                el.kind,
                ElemKind::Hex8 | ElemKind::Hex8I | ElemKind::Hex8R
            ) && fm.face >= 1
            {
                let (ke_f, fe) = heat::hex8_face_heat(&xyz, fm.face, 0.0, fm.h, fm.t_inf)?;
                for a in 0..8 {
                    f[gdofs[a]] += fe[a];
                    for b in 0..8 {
                        let v = ke_f[a * 8 + b];
                        if v.abs() > 0.0 {
                            trips.push((gdofs[a], gdofs[b], v));
                        }
                    }
                }
            }
        }
    }
    for cf in &model.cfluxes {
        let ni = model.node_index(cf.node)?;
        f[ni] += cf.mag;
    }

    let (steady, dt, period) = match model.procedure {
        Procedure::HeatTransfer {
            steady,
            dt,
            period,
        } => (steady, dt, period),
        _ => (true, 0.0, 0.0),
    };

    let mpcs: Vec<crate::constraint::Mpc> = Vec::new();
    let map = DofMap::build(ndof, &prescribed, &mpcs)?;
    let nfree = map.n_ind;
    let mut t_full = map.u0.clone();
    for ic in &model.init {
        if ic.kind == InitKind::Temperature {
            if let Ok(ni) = model.node_index(ic.node) {
                if map.ind_of[ni] >= 0 {
                    t_full[ni] = ic.value;
                }
            }
        }
    }
    let mut solver = "prescribed".to_string();
    let mut residual = 0.0;
    let mut iters = 1usize;
    if nfree > 0 {
        if !steady && dt > 0.0 && period > 0.0 && c_diag.iter().any(|v| *v > 0.0) {
            let nsteps = ((period / dt).round() as usize).clamp(1, 20_000);
            let h = period / nsteps as f64;
            let mut kdt = trips.clone();
            for i in 0..ndof {
                if c_diag[i].abs() > 0.0 {
                    kdt.push((i, i, c_diag[i] / h));
                }
            }
            for _ in 0..nsteps {
                let mut rhs_full = f.clone();
                for i in 0..ndof {
                    rhs_full[i] += c_diag[i] / h * t_full[i];
                }
                let (ff, rhs) = map.reduce(&kdt, &rhs_full);
                let solved = solve_kff(nfree, ff, &rhs)?;
                solver = solved.name;
                residual = solved.residual;
                t_full = map.reconstruct(&solved.x);
            }
            solver = format!("backward-Euler ({solver}, {nsteps} steps)");
            iters = nsteps;
        } else {
            let (ff, rhs) = map.reduce(&trips, &f);
            let solved = solve_kff(nfree, ff, &rhs)?;
            solver = solved.name;
            residual = solved.residual;
            t_full = map.reconstruct(&solved.x);
        }
    }

    let mut kt = vec![0.0; ndof];
    for &(i, j, v) in &trips {
        kt[i] += v * t_full[j];
    }
    let mut rf_full = vec![0.0; ndof];
    for i in 0..ndof {
        rf_full[i] = kt[i] - f[i];
    }

    let mut u = vec![[0.0; 3]; nnode];
    let mut rf = vec![[0.0; 3]; nnode];
    for i in 0..nnode {
        u[i][0] = t_full[i];
        rf[i][0] = rf_full[i];
    }
    let ur = vec![[0.0; 3]; nnode];
    let rm = vec![[0.0; 3]; nnode];
    let stress = vec![[0.0; 6]; nnode];
    let strain = vec![[0.0; 6]; nnode];
    let vm = vec![0.0; nnode];
    let stress_gp = Vec::new();
    let frd_s = frd::write_frd(&model, &u, &stress, &rf, &strain);
    let dat_s = dat::write_dat(&model, &u, &stress_gp, &rf);
    let procedure = model.procedure.name().to_string();
    Ok(SolveOutput {
        model,
        u,
        ur,
        rf,
        rm,
        stress,
        strain,
        von_mises: vm,
        stress_gp,
        frd: frd_s,
        dat: dat_s,
        nfree,
        ndof,
        solver,
        iters,
        residual,
        time_ms: now_ms() - t0,
        procedure,
        frequencies: Vec::new(),
        buckles: Vec::new(),
    })
}

#[cfg(target_arch = "wasm32")]
fn now_ms() -> f64 {
    0.0
}

#[cfg(not(target_arch = "wasm32"))]
fn now_ms() -> f64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs_f64() * 1000.0)
        .unwrap_or(0.0)
}

fn yield_from_curve(curve: &[(f64, f64)], peeq: f64) -> (f64, f64) {
    if curve.is_empty() {
        return (f64::MAX, 0.0);
    }
    if peeq <= curve[0].0 {
        return (curve[0].1, 0.0);
    }
    for w in curve.windows(2) {
        let (p0, s0) = w[0];
        let (p1, s1) = w[1];
        if peeq <= p1 {
            let t = if (p1 - p0).abs() < 1e-16 {
                0.0
            } else {
                (peeq - p0) / (p1 - p0)
            };
            let sy = s0 + t * (s1 - s0);
            let h = if (p1 - p0).abs() < 1e-16 {
                0.0
            } else {
                (s1 - s0) / (p1 - p0)
            };
            return (sy, h.max(0.0));
        }
    }
    let last = curve[curve.len() - 1];
    (last.1, 0.0)
}

fn truss_1d_stress(e: f64, eps: f64, curve: Option<&[(f64, f64)]>) -> (f64, f64) {
    // returns (sigma, tangent E_t)
    let Some(c) = curve else {
        return (e * eps, e);
    };
    if c.is_empty() {
        return (e * eps, e);
    }
    let sign = if eps >= 0.0 { 1.0 } else { -1.0 };
    let aeps = eps.abs();
    // ε = σ/E + peeq, σ = sy(peeq)
    let mut pe = 0.0;
    let (sy0, _) = yield_from_curve(c, 0.0);
    if aeps <= sy0 / e {
        return (sign * e * aeps, e);
    }
    for _ in 0..40 {
        let (sy, h) = yield_from_curve(c, pe);
        let eps_of = sy / e + pe;
        let r = aeps - eps_of;
        if r.abs() < 1e-14 {
            return (sign * sy, (e * h) / (e + h).max(1e-30));
        }
        let den = h / e + 1.0;
        pe += r / den.max(1e-30);
        pe = pe.max(0.0);
    }
    let (sy, h) = yield_from_curve(c, pe);
    (sign * sy, (e * h) / (e + h).max(1e-30))
}

fn solve_truss_newton(model: Model, t0: f64, nlgeom: bool) -> Result<SolveOutput> {
    let ndn = 3;
    let nnode = model.node_ids.len();
    let ndof = ndn * nnode;
    let mut prescribed: HashMap<usize, f64> = HashMap::new();
    for bc in &model.bcs {
        if bc.dof >= ndn {
            continue;
        }
        let ni = model.node_index(bc.node)?;
        prescribed.insert(dof_of(ndn, ni, bc.dof), bc.value);
    }
    let mut f_ext = vec![0.0; ndof];
    for c in &model.cloads {
        if c.dof >= ndn {
            continue;
        }
        let ni = model.node_index(c.node)?;
        f_ext[dof_of(ndn, ni, c.dof)] += c.mag;
    }
    let mpcs = constraint::build_all_mpcs(&model, ndn)?;
    let map = DofMap::build(ndof, &prescribed, &mpcs)?;
    let nfree = map.n_ind;
    if nfree == 0 {
        return err("Truss-Newton: keine freien DOF.");
    }
    let mut u_full = map.u0.clone();
    let mut solver = String::new();
    let mut residual = 0.0;
    let mut iters = 0usize;
    for it in 0..25 {
        iters = it + 1;
        let mut trips: Vec<(usize, usize, f64)> = Vec::new();
        let mut f_int = vec![0.0; ndof];
        for el in &model.elements {
            let xyz0 = elem_xyz(&model, &el.nodes)?;
            let nn = el.kind.nnodes();
            let i1 = if nn == 2 { 1 } else { nn - 1 };
            let mut xyz = xyz0.clone();
            let mut ue = vec![0.0; 3 * nn];
            for a in 0..nn {
                let ni = model.node_index(el.nodes[a])?;
                for d in 0..3 {
                    ue[3 * a + d] = u_full[dof_of(ndn, ni, d)];
                    xyz[a][d] = xyz0[a][d] + ue[3 * a + d];
                }
            }
            let mat = model.material_for(el)?;
            let area = model.thickness_for(el);
            let mut d = [
                xyz[i1][0] - xyz[0][0],
                xyz[i1][1] - xyz[0][1],
                xyz[i1][2] - xyz[0][2],
            ];
            let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt().max(1e-18);
            d[0] /= len;
            d[1] /= len;
            d[2] /= len;
            let mut d0 = [
                xyz0[i1][0] - xyz0[0][0],
                xyz0[i1][1] - xyz0[0][1],
                xyz0[i1][2] - xyz0[0][2],
            ];
            let l0 = (d0[0] * d0[0] + d0[1] * d0[1] + d0[2] * d0[2]).sqrt().max(1e-18);
            d0[0] /= l0;
            d0[1] /= l0;
            d0[2] /= l0;
            let eps = (len - l0) / l0;
            let curve = model.plastic_for(el);
            let (sig, et) = truss_1d_stress(mat.e, eps, curve);
            let nforce = sig * area;
            let kax = et * area / l0;
            let nd = 3 * nn;
            let mut ke = vec![0.0; nd * nd];
            let geom = if nlgeom { nforce / len } else { 0.0 };
            for a in [0usize, i1] {
                for b in [0usize, i1] {
                    let s_ax = if a == b { kax } else { -kax };
                    let s_g = if a == b { geom } else { -geom };
                    for i in 0..3 {
                        for j in 0..3 {
                            let ax = d[i] * d[j];
                            let pr = if i == j { 1.0 } else { 0.0 } - d[i] * d[j];
                            ke[(3 * a + i) * nd + (3 * b + j)] += s_ax * ax + s_g * pr;
                        }
                    }
                }
            }
            let mut fe = vec![0.0; nd];
            for k in 0..3 {
                fe[k] -= nforce * d[k];
                fe[3 * i1 + k] += nforce * d[k];
            }
            if mat.alpha.abs() > 0.0 {
                let t0n = 0.5
                    * (model.temperature_at(el.nodes[0]) + model.temperature_at(el.nodes[i1]));
                let dth = t0n - mat.tref;
                let nth = mat.e * area * mat.alpha * dth;
                for k in 0..3 {
                    fe[k] -= nth * d0[k];
                    fe[3 * i1 + k] += nth * d0[k];
                }
            }
            let mut gdofs = Vec::new();
            for a in 0..nn {
                let ni = model.node_index(el.nodes[a])?;
                for d in 0..3 {
                    gdofs.push(dof_of(ndn, ni, d));
                }
            }
            for i in 0..nd {
                f_int[gdofs[i]] += fe[i];
                for j in 0..nd {
                    let v = ke[i * nd + j];
                    if v.abs() > 0.0 {
                        trips.push((gdofs[i], gdofs[j], v));
                    }
                }
            }
        }
        let mut r = vec![0.0; ndof];
        for i in 0..ndof {
            r[i] = f_ext[i] - f_int[i];
        }
        let (ff, rhs) = map.reduce_inc(&trips, &r);
        let solved = solve_kff(nfree, ff, &rhs)?;
        solver = solved.name;
        residual = rhs.iter().map(|v| v * v).sum::<f64>().sqrt();
        let du_full = {
            let mut t = vec![0.0; ndof];
            // reconstruct increment: T * du (ignore u0)
            for i in 0..ndof {
                for &(j, c) in &map.t_row[i] {
                    t[i] += c * solved.x[j];
                }
            }
            t
        };
        for i in 0..ndof {
            u_full[i] += du_full[i];
        }
        let dun = solved.x.iter().map(|v| v * v).sum::<f64>().sqrt();
        if residual < 1e-8 * (1.0 + f_ext.iter().map(|v| v * v).sum::<f64>().sqrt()) || dun < 1e-12
        {
            break;
        }
        if it == 24 {
            return err(format!(
                "Newton konvergierte nicht (r={residual:.3e} nach {iters} Iterationen)."
            ));
        }
    }
    let mut u = vec![[0.0; 3]; nnode];
    let ur = vec![[0.0; 3]; nnode];
    let mut rf = vec![[0.0; 3]; nnode];
    let rm = vec![[0.0; 3]; nnode];
    // reactions from last residual ~ f_int - f_ext at supports
    let mut f_int = vec![0.0; ndof];
    for el in &model.elements {
        let xyz0 = elem_xyz(&model, &el.nodes)?;
        let nn = el.kind.nnodes();
        let i1 = if nn == 2 { 1 } else { nn - 1 };
        let mut xyz = xyz0.clone();
        for a in 0..nn {
            let ni = model.node_index(el.nodes[a])?;
            for d in 0..3 {
                xyz[a][d] = xyz0[a][d] + u_full[dof_of(ndn, ni, d)];
            }
        }
        let mat = model.material_for(el)?;
        let area = model.thickness_for(el);
        let mut d = [
            xyz[i1][0] - xyz[0][0],
            xyz[i1][1] - xyz[0][1],
            xyz[i1][2] - xyz[0][2],
        ];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt().max(1e-18);
        d[0] /= len;
        d[1] /= len;
        d[2] /= len;
        let l0 = {
            let dx = xyz0[i1][0] - xyz0[0][0];
            let dy = xyz0[i1][1] - xyz0[0][1];
            let dz = xyz0[i1][2] - xyz0[0][2];
            (dx * dx + dy * dy + dz * dz).sqrt().max(1e-18)
        };
        let eps = (len - l0) / l0;
        let (sig, _) = truss_1d_stress(mat.e, eps, model.plastic_for(el));
        let nforce = sig * area;
        for a in 0..nn {
            let ni = model.node_index(el.nodes[a])?;
            let s = if a == 0 { -1.0 } else if a == i1 { 1.0 } else { 0.0 };
            for k in 0..3 {
                f_int[dof_of(ndn, ni, k)] += s * nforce * d[k];
            }
        }
    }
    for ni in 0..nnode {
        for d in 0..3 {
            u[ni][d] = u_full[dof_of(ndn, ni, d)];
            rf[ni][d] = f_int[dof_of(ndn, ni, d)] - f_ext[dof_of(ndn, ni, d)];
        }
    }
    let mut stress = vec![[0.0; 6]; nnode];
    let mut strain = vec![[0.0; 6]; nnode];
    let mut vm = vec![0.0; nnode];
    let mut stress_gp = Vec::new();
    for el in &model.elements {
        let xyz0 = elem_xyz(&model, &el.nodes)?;
        let nn = el.kind.nnodes();
        let i1 = if nn == 2 { 1 } else { nn - 1 };
        let mut xyz = xyz0.clone();
        for a in 0..nn {
            let ni = model.node_index(el.nodes[a])?;
            for d in 0..3 {
                xyz[a][d] += u[ni][d];
            }
        }
        let mat = model.material_for(el)?;
        let l0 = {
            let dx = xyz0[i1][0] - xyz0[0][0];
            let dy = xyz0[i1][1] - xyz0[0][1];
            let dz = xyz0[i1][2] - xyz0[0][2];
            (dx * dx + dy * dy + dz * dz).sqrt().max(1e-18)
        };
        let len = {
            let dx = xyz[i1][0] - xyz[0][0];
            let dy = xyz[i1][1] - xyz[0][1];
            let dz = xyz[i1][2] - xyz[0][2];
            (dx * dx + dy * dy + dz * dz).sqrt().max(1e-18)
        };
        let eps = (len - l0) / l0;
        let (sig, _) = truss_1d_stress(mat.e, eps, model.plastic_for(el));
        let mut d = [
            xyz[i1][0] - xyz[0][0],
            xyz[i1][1] - xyz[0][1],
            xyz[i1][2] - xyz[0][2],
        ];
        d[0] /= len;
        d[1] /= len;
        d[2] /= len;
        let s = [
            sig * d[0] * d[0],
            sig * d[1] * d[1],
            sig * d[2] * d[2],
            sig * d[0] * d[1],
            sig * d[1] * d[2],
            sig * d[2] * d[0],
        ];
        for a in 0..nn {
            let ni = model.node_index(el.nodes[a])?;
            stress[ni] = s;
            strain[ni][0] = eps * d[0] * d[0];
            vm[ni] = sig.abs();
        }
        stress_gp.push((el.id, 1usize, s));
    }
    let frd_s = frd::write_frd(&model, &u, &stress, &rf, &strain);
    let dat_s = dat::write_dat(&model, &u, &stress_gp, &rf);
    let procedure = model.procedure.name().to_string();
    Ok(SolveOutput {
        model,
        u,
        ur,
        rf,
        rm,
        stress,
        strain,
        von_mises: vm,
        stress_gp,
        frd: frd_s,
        dat: dat_s,
        nfree,
        ndof,
        solver,
        iters,
        residual,
        time_ms: now_ms() - t0,
        procedure,
        frequencies: vec![],
        buckles: vec![],
    })
}
