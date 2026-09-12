use std::collections::HashMap;

use crate::beam;
use crate::dat;
use crate::elem::{
    element_ke, element_nodal_stress, hex8_body_force, hex8_face_pressure, quad4_body_force,
    quad4_edge_pressure, tet4_body_force, von_mises,
};
use crate::error::{err, Result};
use crate::frd;
use crate::linalg::{chol_solve, csr_from_triplets, pcg};
use crate::model::{Dload, ElemKind, Model};
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

    let mut is_free = vec![true; ndof];
    for (&d, _) in &prescribed {
        if d < ndof {
            is_free[d] = false;
        }
    }
    let mut free_of = vec![-1isize; ndof];
    let mut nfree = 0usize;
    for d in 0..ndof {
        if is_free[d] {
            free_of[d] = nfree as isize;
            nfree += 1;
        }
    }
    if nfree == 0 {
        return err("Alle Freiheitsgrade sind gelagert — nichts zu lösen.");
    }

    let mut f_full = vec![0.0; ndof];
    for c in &model.cloads {
        if c.dof >= ndn {
            continue;
        }
        let ni = model.node_index(c.node)?;
        f_full[dof_of(ndn, ni, c.dof)] += c.mag;
    }

    let mut trips: Vec<(usize, usize, f64)> = Vec::new();
    for el in &model.elements {
        let xyz = elem_xyz(&model, &el.nodes)?;
        let mat = model.material_for(el)?;
        let th = model.thickness_for(el);
        let sec = if el.kind.is_beam() {
            Some(model.beam_section_for(el)?)
        } else {
            None
        };
        let kef = element_ke(el.kind, &xyz, mat.e, mat.nu, th, sec.as_ref())?;
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
    }

    // Reduce to free system: Kff u_f = f_f - Kfp u_p
    let mut rhs = vec![0.0; nfree];
    for d in 0..ndof {
        if is_free[d] {
            rhs[free_of[d] as usize] += f_full[d];
        }
    }
    let mut ff_trips: Vec<(usize, usize, f64)> = Vec::new();
    for (i, j, v) in &trips {
        let fi = free_of[*i];
        let fj = free_of[*j];
        if fi >= 0 && fj >= 0 {
            ff_trips.push((fi as usize, fj as usize, *v));
        } else if fi >= 0 && fj < 0 {
            let up = prescribed.get(j).copied().unwrap_or(0.0);
            rhs[fi as usize] -= *v * up;
        }
    }

    const DENSE_LIMIT: usize = 900;
    let (u_free, solver, iters, residual) = if nfree <= DENSE_LIMIT {
        let csr = csr_from_triplets(nfree, ff_trips);
        let mut dense = csr.to_dense();
        let x = chol_solve(&mut dense, nfree, &rhs)?;
        let mut r = vec![0.0; nfree];
        csr.matvec(&x, &mut r);
        let mut res = 0.0;
        for i in 0..nfree {
            let d = r[i] - rhs[i];
            res += d * d;
        }
        (x, "Cholesky".to_string(), 1usize, res.sqrt())
    } else {
        let csr = csr_from_triplets(nfree, ff_trips);
        let (x, info) = pcg(&csr, &rhs, 1e-8, (4 * nfree).max(200))?;
        (x, "PCG".to_string(), info.iters, info.residual)
    };

    let mut u_full = vec![0.0; ndof];
    for d in 0..ndof {
        if is_free[d] {
            u_full[d] = u_free[free_of[d] as usize];
        } else if let Some(&v) = prescribed.get(&d) {
            u_full[d] = v;
        }
    }

    // Reactions: R = K u - F_applied (nonzero on supports)
    let mut ku = vec![0.0; ndof];
    for (i, j, v) in &trips {
        ku[*i] += *v * u_full[*j];
    }
    let mut rf_full = vec![0.0; ndof];
    for d in 0..ndof {
        rf_full[d] = ku[d] - f_full[d];
    }

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

    let frd_s = frd::write_frd(&model, &u, &stress, &rf, &strain);
    let dat_s = dat::write_dat(&model, &u, &stress_gp, &rf);
    let dt = now_ms() - t0;

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
        ElemKind::Hex8 => {
            let mut p = [[0.0; 3]; 8];
            for i in 0..8 {
                p[i] = xyz[i];
            }
            let fe = hex8_face_pressure(&p, face, mag)?;
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
        ElemKind::Hex8 => {
            let mut p = [[0.0; 3]; 8];
            for i in 0..8 {
                p[i] = xyz[i];
            }
            let fe = hex8_body_force(&p, bx, by, bz)?;
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
