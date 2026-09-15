use std::collections::HashMap;
use std::io::Write;

use crate::beam;
use crate::constraint::{self, DofMap};
use crate::contact;
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
use crate::material::truss_1d_stress;
use crate::model::{Dload, ElemKind, FluxKind, InitKind, Material, Model, Procedure, RiksCtrl};
use crate::nlgeom;
use crate::plastic;
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
    pub nsteps: usize,
    pub peeq: Vec<f64>,
    pub lambda: f64,
    pub ninc: usize,
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

/// Limit each node's displacement increment so a mechanism in one DOF cannot
/// starve the pretension dummy (or invert a contact pair) in the same step.
fn cap_nodal_du(du: &mut [f64], ndn: usize, nnode: usize, cap: f64) {
    if cap <= 0.0 || !cap.is_finite() {
        return;
    }
    let dim = 3.min(ndn);
    for ni in 0..nnode {
        let base = ndn * ni;
        if base >= du.len() {
            break;
        }
        let mut nrm = 0.0_f64;
        for d in 0..dim {
            let i = base + d;
            if i >= du.len() {
                break;
            }
            if !du[i].is_finite() {
                du[i] = 0.0;
                nrm = cap * 2.0;
            } else {
                nrm = nrm.max(du[i].abs());
            }
        }
        if nrm > cap {
            let s = cap / nrm;
            for d in 0..dim {
                let i = base + d;
                if i < du.len() {
                    du[i] *= s;
                }
            }
        }
    }
}

pub fn solve(model: Model) -> Result<SolveOutput> {
    let t0 = now_ms();
    if model.steps.len() > 1 {
        return solve_sequence(model, t0);
    }
    solve_one(model, t0)
}

fn solve_sequence(model: Model, t0: f64) -> Result<SolveOutput> {
    let steps = model.steps.clone();
    let all_c = model.cloads.clone();
    let all_d = model.dloads.clone();
    let all_b = model.bcs.clone();
    let mut last: Option<SolveOutput> = None;
    let mut u_prev: Vec<[f64; 3]> = Vec::new();
    let mut f_prev: Vec<f64> = Vec::new();
    for st in &steps {
        let mut m = model.clone();
        m.procedure = st.procedure.clone();
        let c0 = st.cload_from.min(all_c.len());
        let c1 = st.n_cload.min(all_c.len()).max(c0);
        let d0 = st.dload_from.min(all_d.len());
        let d1 = st.n_dload.min(all_d.len()).max(d0);
        let b0 = st.bc_from.min(all_b.len());
        let b1 = st.n_bc.min(all_b.len()).max(b0);
        m.cloads = all_c[c0..c1].to_vec();
        m.dloads = all_d[d0..d1].to_vec();
        m.bcs = all_b[b0..b1].to_vec();
        m.steps.clear();
        m.u_start = u_prev.clone();
        m.f_start = f_prev.clone();
        let mut out = solve_one(m, t0)?;
        out.nsteps = steps.len();
        u_prev = out.u.clone();
        let ndn = out.model.ndof_node();
        let ndof = ndn * out.model.node_ids.len();
        f_prev = assemble_fext(&out.model, ndn, ndof, out.model.static_period)?;
        last = Some(out);
    }
    last.ok_or_else(|| crate::error::FemError("Keine Schritte.".into()))
}

fn solve_one(mut model: Model, t0: f64) -> Result<SolveOutput> {
    remap_rot_node_bcs(&mut model);
    if matches!(model.procedure, Procedure::HeatTransfer { .. }) {
        return solve_heat(model, t0);
    }
    let nlgeom = matches!(
        model.procedure,
        Procedure::Static { nlgeom: true, .. }
    );
    let riks = matches!(model.procedure, Procedure::Static { riks: true, .. });
    let truss2_only = !model.elements.is_empty()
        && model.elements.iter().all(|e| e.kind == ElemKind::Truss2);
    let continuum_only = !model.elements.is_empty()
        && model.elements.iter().all(|e| nlgeom::is_nl_continuum(e.kind));
    let riks_ok = !model.elements.is_empty()
        && model
            .elements
            .iter()
            .all(|e| e.kind.is_truss() || nlgeom::is_nl_continuum(e.kind));
    if model.has_contact() {
        if !matches!(model.procedure, Procedure::Static { .. }) {
            model.warn(
                "*CONTACT PAIR mit nicht-statischer Prozedur: lineare Kontaktlösung.",
            );
        }
        if model.has_plastic() || riks {
            model.warn(
                "*CONTACT PAIR mit NLGEOM: Kontakt auf deformierter Geometrie im Newton.",
            );
            return solve_continuum_newton(model, t0);
        }
        if nlgeom {
            model.warn(
                "*CONTACT PAIR mit NLGEOM: Kontakt-Newton auf linearisierter Steifigkeit.",
            );
        }
        let mut fallback = model.clone();
        fallback.contact_pairs.clear();
        match solve_contact(model, t0) {
            Ok(o) => return Ok(o),
            Err(e) => {
                if nlgeom {
                    return Err(e);
                }
                fallback.warn(format!(
                    "Kontakt nicht konvergiert ({e}) — linear ohne Kontakt."
                ));
                return solve_linear(fallback, t0);
            }
        }
    }
    if riks {
        if !riks_ok {
            return err("RIKS ist für T3D2 und Kontinuum (C3D*) implementiert; gemischte Netze nicht.");
        }
        return solve_riks(model, t0);
    }
    if nlgeom || model.has_plastic() {
        if truss2_only {
            return solve_truss_newton(model, t0, nlgeom);
        }
        if continuum_only && model.has_plastic() && !nlgeom {
            return solve_continuum_plastic(model, t0);
        }
        if nlgeom && can_nlgeom_newton(&model) {
            return solve_continuum_newton(model, t0);
        }
        if continuum_only && model.has_plastic() {
            return solve_continuum_plastic(model, t0);
        }
        model.warn(
            "NLGEOM/*PLASTIC für diesen Elementmix nicht verfügbar — linear-elastisch gerechnet.",
        );
        return solve_linear(model, t0);
    }
    solve_linear(model, t0)
}

fn can_nlgeom_newton(model: &Model) -> bool {
    !model.elements.is_empty()
        && model.elements.iter().all(|e| {
            nlgeom::is_nl_continuum(e.kind)
                || e.kind.is_truss()
                || e.kind.is_beam()
                || e.kind.is_shell()
                || e.kind.is_membrane()
                || e.kind.is_special()
                || e.kind.is_spring()
        })
}

fn assemble_nl_element(
    model: &Model,
    el: &crate::model::Element,
    ei: usize,
    xyz0: &[[f64; 3]],
    ue: &[f64],
    hist: &[Vec<plastic::GpHist>],
    trial_hist: &mut [Vec<plastic::GpHist>],
) -> Result<nlgeom::NlElem> {
    if el.kind.is_beam() {
        let mat = model.material_for(el)?;
        let sec = model.beam_section_for(el)?;
        let (ke, fe, cauchy) = beam::stiffness_nl(el.kind, xyz0, ue, mat.e, mat.nu, &sec)?;
        return Ok(nlgeom::NlElem {
            ke,
            fe,
            vol: 1.0,
            cauchy,
            gl: [0.0; 6],
            peeq: 0.0,
        });
    }
    if el.kind.is_shell() {
        let mat = model.material_for(el)?;
        let th = model.thickness_for(el);
        let (ke, fe, cauchy) = shell::stiffness_nl(el.kind, xyz0, ue, mat.e, mat.nu, th)?;
        return Ok(nlgeom::NlElem {
            ke,
            fe,
            vol: 1.0,
            cauchy,
            gl: [0.0; 6],
            peeq: 0.0,
        });
    }
    if el.kind.is_truss() {
        return truss_nl_elem(model, el, xyz0, ue);
    }
    if el.kind.is_spring() || el.kind.is_membrane() {
        let mat = if el.kind.needs_material() {
            model.material_for(el)?
        } else {
            Material::default()
        };
        let th = if el.kind.is_spring() {
            model.spring_k_for(el)?
        } else {
            model.thickness_for(el)
        };
        let kef = element_ke(el.kind, xyz0, mat.e, mat.nu, th, None)?;
        let n = kef.ndof;
        let mut fe = vec![0.0; n];
        for i in 0..n {
            for j in 0..n.min(ue.len()) {
                fe[i] += kef.ke[i * n + j] * ue[j];
            }
        }
        return Ok(nlgeom::NlElem {
            ke: kef.ke,
            fe,
            vol: 1.0,
            cauchy: [0.0; 6],
            gl: [0.0; 6],
            peeq: 0.0,
        });
    }
    let mat = model.material_for(el)?;
    let th = model.thickness_for(el);
    if el.kind.is_continuum3d() {
        if let Some(curve) = model.plastic_for(el) {
            let empty: Vec<plastic::GpHist> = Vec::new();
            let h = hist.get(ei).map(|v| v.as_slice()).unwrap_or(&empty);
            let (n, hnew) =
                plastic::continuum_plastic_nl(el.kind, xyz0, ue, mat.e, mat.nu, curve, h)?;
            if ei < trial_hist.len() {
                trial_hist[ei] = hnew;
            }
            return Ok(n);
        }
    }
    nlgeom::continuum_nl(el.kind, xyz0, ue, &mat, th)
}

fn truss_nl_elem(
    model: &Model,
    el: &crate::model::Element,
    xyz0: &[[f64; 3]],
    ue: &[f64],
) -> Result<nlgeom::NlElem> {
    let nn = el.kind.nnodes();
    let i1 = if nn == 2 { 1 } else { nn - 1 };
    let mut xyz = xyz0.to_vec();
    for a in 0..nn {
        for d in 0..3 {
            xyz[a][d] = xyz0[a][d] + ue.get(3 * a + d).copied().unwrap_or(0.0);
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
    let _ = d0;
    let eps = (len - l0) / l0;
    let (sig, et) = truss_1d_stress(mat.e, eps, model.plastic_for(el));
    let nforce = sig * area;
    let kax = et * area / l0;
    let nd = 3 * nn;
    let mut ke = vec![0.0; nd * nd];
    let geom = nforce / len;
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
    Ok(nlgeom::NlElem {
        ke,
        fe,
        vol: l0 * area,
        cauchy: [
            sig * d[0] * d[0],
            sig * d[1] * d[1],
            sig * d[2] * d[2],
            sig * d[0] * d[1],
            sig * d[1] * d[2],
            sig * d[2] * d[0],
        ],
        gl: [eps, 0.0, 0.0, 0.0, 0.0, 0.0],
        peeq: 0.0,
    })
}

fn solve_linear(model: Model, t0: f64) -> Result<SolveOutput> {
    let ndn = model.ndof_node();
    let nnode = model.node_ids.len();
    let ndof = ndn * nnode;
    if ndof == 0 {
        return err("Modell ohne Freiheitsgrade.");
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
    pin_unused_dofs(&model, ndn, nnode, &mut prescribed)?;

    let mut f_full = vec![0.0; ndof];
    let mut m_full = vec![0.0; ndof];

    let mut trips: Vec<(usize, usize, f64)> = Vec::new();
    let mut c_trips: Vec<(usize, usize, f64)> = Vec::new();
    for el in &model.elements {
        let xyz = elem_xyz(&model, &el.nodes)?;
        if el.kind.is_special() {
            scatter_special(
                &model,
                el,
                &xyz,
                ndn,
                &mut m_full,
                &mut trips,
                &mut c_trips,
            )?;
            apply_point_grav(&model, el, ndn, &mut f_full)?;
            continue;
        }
        let mat = if el.kind.needs_material() {
            model.material_for(el)?
        } else {
            Material::default()
        };
        let th = if el.kind.is_spring() {
            model.spring_k_for(el)?
        } else if el.kind.is_gap() {
            model.gap_for(el)?.k
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
                Dload::Centrif {
                    omega2,
                    p1,
                    p2,
                    elems,
                } if centrif_applies(elems, el.id) => {
                    let a = centrif_accel(*omega2, *p1, *p2, &xyz);
                    let bx = mat.density * a[0];
                    let by = mat.density * a[1];
                    let bz = mat.density * a[2];
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

    let t_cload = if matches!(model.procedure, Procedure::Dynamic { .. }) {
        0.0
    } else {
        model.static_period.max(0.0)
    };
    add_cloads(&model, ndn, t_cload, &mut f_full)?;
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
            &c_trips,
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
        if !el.kind.needs_material() {
            continue;
        }
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
                } else if matches!(el.kind, ElemKind::Hex20 | ElemKind::Hex20R) {
                    let mut mean = [0.0; 6];
                    if let Some((_, _, s)) = stress_gp.iter().find(|(id, _, _)| *id == el.id) {
                        mean = *s;
                    }
                    crate::eigen::hex20_kg(&xyz, &mean, el.kind.reduced_int())?
                } else if el.kind == ElemKind::Tet4 {
                    let mut mean = [0.0; 6];
                    if let Some((_, _, s)) = stress_gp.iter().find(|(id, _, _)| *id == el.id) {
                        mean = *s;
                    }
                    crate::eigen::tet4_kg(&xyz, &mean)?
                } else if el.kind == ElemKind::Tet10 || el.kind == ElemKind::Tet10T {
                    let mut mean = [0.0; 6];
                    if let Some((_, _, s)) = stress_gp.iter().find(|(id, _, _)| *id == el.id) {
                        mean = *s;
                    }
                    crate::eigen::tet10_kg(&xyz, &mean)?
                } else if el.kind == ElemKind::Wedge6 {
                    let mut mean = [0.0; 6];
                    if let Some((_, _, s)) = stress_gp.iter().find(|(id, _, _)| *id == el.id) {
                        mean = *s;
                    }
                    extra::wedge6_kg(&xyz, &mean)?
                } else if el.kind == ElemKind::Wedge15 {
                    let mut mean = [0.0; 6];
                    if let Some((_, _, s)) = stress_gp.iter().find(|(id, _, _)| *id == el.id) {
                        mean = *s;
                    }
                    extra::wedge15_kg(&xyz, &mean)?
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

    let frd_s = frd::write_frd(&model, &u, &stress, &rf, &strain, &[]);
    let dat_s = dat::write_dat(&model, &u, &stress_gp, &rf);
    let dt = now_ms() - t0;
    let procedure = model.procedure.name().to_string();
    let nsteps = model.steps.len().max(1);

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
        nsteps,
        lambda: 1.0,
        ninc: 1,
        peeq: Vec::new(),
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
        ElemKind::Wedge15 => {
            let fe = extra::wedge15_face_pressure(xyz, face, mag)?;
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
        ElemKind::Cax4 | ElemKind::Cax4R | ElemKind::Cax8 | ElemKind::Cax8R | ElemKind::Cax3 | ElemKind::Cax6 => {
            let fe = crate::axisym::cax_edge_pressure(xyz, kind.nnodes(), face, mag)?;
            scatter_fe(&fe, gdofs, 2, local_dim, f_full);
        }
        ElemKind::Mem3 | ElemKind::Mem4 | ElemKind::Mem4R | ElemKind::Mem6 | ElemKind::Mem8 => {
            let fe = shell::membrane_pressure(kind, xyz, mag)?;
            scatter_fe(&fe, gdofs, 3, local_dim, f_full);
        }
        k if k.is_shell() => {
            let fe = shell::pressure_force(k, xyz, mag)?;
            scatter_fe(&fe, gdofs, 6, local_dim, f_full);
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
        ElemKind::Wedge15 => {
            let fe = extra::wedge15_body_force(xyz, bx, by, bz)?;
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
        ElemKind::Tet10 | ElemKind::Tet10T => {
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
        ElemKind::Cax4 | ElemKind::Cax4R => {
            let fe = crate::axisym::cax4_body_force(xyz, bx, by, kind.reduced_int())?;
            scatter_fe(&fe, gdofs, 2, local_dim, f_full);
        }
        ElemKind::Cax8 | ElemKind::Cax8R => {
            let fe = crate::axisym::cax8_body_force(xyz, bx, by, kind.reduced_int())?;
            scatter_fe(&fe, gdofs, 2, local_dim, f_full);
        }
        ElemKind::Cax3 => {
            let fe = crate::axisym::cax3_body_force(xyz, bx, by)?;
            scatter_fe(&fe, gdofs, 2, local_dim, f_full);
        }
        ElemKind::Cax6 => {
            let fe = crate::axisym::cax6_body_force(xyz, bx, by)?;
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

/// Centrifugal acceleration at the element centroid: ω² r_⊥.
fn centrif_accel(omega2: f64, p1: [f64; 3], p2: [f64; 3], xyz: &[[f64; 3]]) -> [f64; 3] {
    if xyz.is_empty() {
        return [0.0, 0.0, 0.0];
    }
    let n = xyz.len() as f64;
    let mut c = [0.0; 3];
    for p in xyz {
        c[0] += p[0];
        c[1] += p[1];
        c[2] += p[2];
    }
    c[0] /= n;
    c[1] /= n;
    c[2] /= n;
    let mut axis = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    let al = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if al < 1e-18 {
        axis = [0.0, 0.0, 1.0];
    } else {
        axis[0] /= al;
        axis[1] /= al;
        axis[2] /= al;
    }
    let r = [c[0] - p1[0], c[1] - p1[1], c[2] - p1[2]];
    let proj = r[0] * axis[0] + r[1] * axis[1] + r[2] * axis[2];
    let rp = [r[0] - proj * axis[0], r[1] - proj * axis[1], r[2] - proj * axis[2]];
    [omega2 * rp[0], omega2 * rp[1], omega2 * rp[2]]
}

fn centrif_applies(elems: &[i32], id: i32) -> bool {
    elems.is_empty() || elems.contains(&id)
}

/// CalculiX *RIGID BODY, ROT NODE=n: dofs 1–3 of the rot node are θ of the ref node.
fn remap_rot_node_bcs(model: &mut Model) {
    let maps: Vec<(i32, i32)> = model
        .rigid_bodies
        .iter()
        .filter_map(|rb| rb.rot_node.map(|r| (r, rb.ref_node)))
        .collect();
    if maps.is_empty() {
        return;
    }
    for bc in &mut model.bcs {
        for &(rot, refn) in &maps {
            if bc.node == rot && bc.dof < 3 {
                bc.node = refn;
                bc.dof += 3;
                break;
            }
        }
    }
}

/// Pin unused rotational DOFs on continuum nodes, and uz on planar 2D/truss models.
fn pin_unused_dofs(
    model: &Model,
    ndn: usize,
    nnode: usize,
    prescribed: &mut HashMap<usize, f64>,
) -> Result<()> {
    if ndn == 6 {
        let mut struct_node = vec![false; nnode];
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
        for ni in 0..nnode {
            if !struct_node[ni] {
                for r in 3..6 {
                    prescribed.entry(dof_of(ndn, ni, r)).or_insert(0.0);
                }
            }
        }
    }
    // Isolated nodes (dummy REF/ROT nodes, unused mesh ids) get all DOFs pinned
    // unless they are a rigid/coupling reference.
    {
        let mut attached = vec![false; nnode];
        for el in &model.elements {
            for &id in &el.nodes {
                if let Ok(i) = model.node_index(id) {
                    attached[i] = true;
                }
            }
        }
        let mut keep = vec![false; nnode];
        for rb in &model.rigid_bodies {
            if let Ok(i) = model.node_index(rb.ref_node) {
                keep[i] = true;
            }
        }
        for c in &model.couplings {
            if let Ok(i) = model.node_index(c.ref_node) {
                keep[i] = true;
            }
        }
        for p in &model.pretensions {
            if let Ok(i) = model.node_index(p.dummy) {
                keep[i] = true;
                // Dummy dof 1 is the pretension DOF (MPC slave / CLOAD target).
                // Pin the remaining dofs so the isolated node is not singular.
                for d in 1..ndn {
                    prescribed.entry(dof_of(ndn, i, d)).or_insert(0.0);
                }
            }
        }
        for ni in 0..nnode {
            if !attached[ni] && !keep[ni] {
                for d in 0..ndn {
                    prescribed.entry(dof_of(ndn, ni, d)).or_insert(0.0);
                }
            }
        }
    }
    if ndn >= 3 {
        let z0 = model.coords.first().map(|c| c[2]).unwrap_or(0.0);
        let planar = !model.coords.is_empty()
            && model.coords.iter().all(|c| (c[2] - z0).abs() <= 1e-12);
        let has_solid = model.elements.iter().any(|e| e.kind.is_continuum3d());
        if planar && !has_solid && !model.has_beams() && !model.has_shells() {
            for ni in 0..nnode {
                prescribed.entry(dof_of(ndn, ni, 2)).or_insert(0.0);
            }
        }
    }
    Ok(())
}

fn nl_eprint(msg: &str) {
    #[cfg(not(target_arch = "wasm32"))]
    {
        let _ = writeln!(std::io::stderr(), "{msg}");
        let _ = std::io::stderr().flush();
    }
    #[cfg(target_arch = "wasm32")]
    {
        let _ = msg;
    }
}

/// CalculiX `nonlingeo.c` / `checkconvergence.c` increment header.
fn nl_print_increment(
    iinc: usize,
    attempt: usize,
    dt_time: f64,
    prev_time: f64,
    step_time: f64,
    total_time: f64,
) {
    nl_eprint(&format!(" increment {iinc} attempt {attempt} "));
    nl_eprint(&format!(" increment size= {dt_time:e}"));
    nl_eprint(&format!(" sum of previous increments={prev_time:e}"));
    nl_eprint(&format!(" actual step time={step_time:e}"));
    nl_eprint(&format!(" actual total time={total_time:e}"));
}

fn nl_print_iteration(
    model: &Model,
    ndn: usize,
    it: usize,
    r: &[f64],
    f_int: &[f64],
    u_inc: &[f64],
    du: &[f64],
    prescribed: &HashMap<usize, f64>,
) {
    nl_eprint("");
    nl_eprint(&format!(" iteration {it}"));
    nl_eprint("");
    let mut n_force = 0.0;
    let mut sum_force = 0.0;
    for &v in f_int {
        let a = v.abs();
        if a > 1e-20 {
            sum_force += a;
            n_force += 1.0;
        }
    }
    let avg = if n_force > 0.0 {
        sum_force / n_force
    } else {
        0.0
    };
    nl_eprint(&format!(" average force= {avg}"));
    nl_eprint(&format!(" time avg. forc= {avg}"));
    let mut rmax = 0.0;
    let mut r_dof = 0usize;
    for (i, &v) in r.iter().enumerate() {
        if prescribed.contains_key(&i) {
            continue;
        }
        let a = v.abs();
        if a >= rmax {
            rmax = a;
            r_dof = i;
        }
    }
    let (rn, rd) = dof_to_node(model, ndn, r_dof);
    nl_eprint(&format!(
        " largest residual force= {rmax} in node {rn} and dof {rd}"
    ));
    let mut umax = 0.0_f64;
    for &v in u_inc {
        umax = umax.max(v.abs());
    }
    nl_eprint(&format!(" largest increment of disp= {umax:e}"));
    let mut dmax = 0.0;
    let mut d_dof = 0usize;
    for (i, &v) in du.iter().enumerate() {
        let a = v.abs();
        if a >= dmax {
            dmax = a;
            d_dof = i;
        }
    }
    let (dn, dd) = dof_to_node(model, ndn, d_dof);
    nl_eprint(&format!(
        " largest correction to disp= {dmax:e} in node {dn} and dof {dd}"
    ));
}

fn dof_to_node(model: &Model, ndn: usize, dof: usize) -> (i32, usize) {
    if ndn == 0 || model.node_ids.is_empty() {
        return (0, 1);
    }
    let ni = dof / ndn;
    let d = dof % ndn + 1;
    let id = model.node_ids.get(ni).copied().unwrap_or(0);
    (id, d)
}

fn average_abs_force(f: &[f64]) -> f64 {
    let mut n = 0.0;
    let mut s = 0.0;
    for &v in f {
        let a = v.abs();
        if a > 1e-20 {
            s += a;
            n += 1.0;
        }
    }
    if n > 0.0 {
        s / n
    } else {
        0.0
    }
}

/// GRAV on *MASS: F = m · mag · dir̂. Concentrated mass has no continuum density.
fn apply_point_grav(
    model: &Model,
    el: &crate::model::Element,
    ndn: usize,
    f: &mut [f64],
) -> Result<()> {
    if el.kind != ElemKind::Mass {
        return Ok(());
    }
    let m = model.mass_for(el)?;
    let ni = model.node_index(el.nodes[0])?;
    for dl in &model.dloads {
        if let Dload::Grav { mag, dir } = dl {
            let mut ndir = *dir;
            let len = (ndir[0] * ndir[0] + ndir[1] * ndir[1] + ndir[2] * ndir[2]).sqrt();
            if len > 0.0 {
                ndir[0] /= len;
                ndir[1] /= len;
                ndir[2] /= len;
            }
            for d in 0..3.min(ndn) {
                f[dof_of(ndn, ni, d)] += m * *mag * ndir[d];
            }
        }
    }
    Ok(())
}

fn scatter_special(
    model: &Model,
    el: &crate::model::Element,
    xyz: &[[f64; 3]],
    ndn: usize,
    m_full: &mut [f64],
    trips: &mut Vec<(usize, usize, f64)>,
    c_trips: &mut Vec<(usize, usize, f64)>,
) -> Result<()> {
    if el.kind == ElemKind::Mass {
        let m = model.mass_for(el)?;
        let ni = model.node_index(el.nodes[0])?;
        for d in 0..3.min(ndn) {
            m_full[dof_of(ndn, ni, d)] += m;
        }
        return Ok(());
    }
    if el.kind == ElemKind::RotaryI {
        let ijk = model.rotary_for(el)?;
        let ni = model.node_index(el.nodes[0])?;
        if ndn >= 6 {
            m_full[dof_of(ndn, ni, 3)] += ijk[0];
            m_full[dof_of(ndn, ni, 4)] += ijk[1];
            m_full[dof_of(ndn, ni, 5)] += ijk[2];
        }
        return Ok(());
    }
    if el.kind == ElemKind::DashpotA {
        let c = model.dashpot_for(el)?;
        let (ke, _) = extra::spring_stiffness(xyz, c)?;
        let n0 = model.node_index(el.nodes[0])?;
        let n1 = model.node_index(el.nodes[1])?;
        let gd = [
            dof_of(ndn, n0, 0),
            dof_of(ndn, n0, 1),
            dof_of(ndn, n0, 2),
            dof_of(ndn, n1, 0),
            dof_of(ndn, n1, 1),
            dof_of(ndn, n1, 2),
        ];
        for i in 0..6 {
            for j in 0..6 {
                let v = ke[i * 6 + j];
                if v.abs() > 0.0 {
                    c_trips.push((gd[i], gd[j], v));
                }
            }
        }
        return Ok(());
    }
    if el.kind == ElemKind::GapUni {
        let g = model.gap_for(el)?;
        if g.clearance > 1e-14 {
            return Ok(());
        }
        let (ke, _) = extra::spring_stiffness(xyz, g.k)?;
        let n0 = model.node_index(el.nodes[0])?;
        let n1 = model.node_index(el.nodes[1])?;
        let gd = [
            dof_of(ndn, n0, 0),
            dof_of(ndn, n0, 1),
            dof_of(ndn, n0, 2),
            dof_of(ndn, n1, 0),
            dof_of(ndn, n1, 1),
            dof_of(ndn, n1, 2),
        ];
        for i in 0..6 {
            for j in 0..6 {
                let v = ke[i * 6 + j];
                if v.abs() > 0.0 {
                    trips.push((gd[i], gd[j], v));
                }
            }
        }
    }
    Ok(())
}

fn add_cloads(model: &Model, ndn: usize, t: f64, f: &mut [f64]) -> Result<()> {
    let period = model.static_period.max(1e-30);
    let pretension_dummy: std::collections::HashSet<i32> =
        model.pretensions.iter().map(|p| p.dummy).collect();
    for c in &model.cloads {
        if c.dof >= ndn {
            continue;
        }
        let ni = model.node_index(c.node)?;
        let s = if c.amplitude.is_empty() && pretension_dummy.contains(&c.node) {
            // Ramp pretension with step time so DIRECT cutbacks actually
            // reduce the bolt load. Named amplitudes still follow their cards.
            (t / period).clamp(0.0, 1.0)
        } else {
            model.amp_value(&c.amplitude, t)
        };
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
    c_trips: &[(usize, usize, f64)],
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
    for &(i, j, v) in c_trips {
        keff_trips.push((i, j, v * a1));
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
    let mut cv = vec![0.0; ndof];
    for &(i, j, v) in c_trips {
        cv[i] += v * v_full[j];
    }
    let mut a_full = vec![0.0; ndof];
    for i in 0..ndof {
        if map.ind_of[i] < 0 {
            v_full[i] = 0.0;
            continue;
        }
        let rhs = f0[i] - ku[i] - alpha_r * m_full[i] * v_full[i] - beta_r * kv[i] - cv[i];
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
        let mut cpred = vec![0.0; ndof];
        for &(i, j, v) in c_trips {
            cpred[i] += v * pred[j];
        }
        for i in 0..ndof {
            reff[i] += cpred[i];
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
        if !el.kind.needs_material() {
            continue;
        }
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
    let frd_s = frd::write_frd(&model, &u, &stress, &rf, &strain, &[]);
    let dat_s = dat::write_dat(&model, &u, &stress_gp, &rf);
    let procedure = model.procedure.name().to_string();
    let nsteps = model.steps.len().max(1);
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
        nsteps,
        lambda: 1.0,
        ninc: 1,
        peeq: Vec::new(),
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

fn assemble_fext(model: &Model, ndn: usize, ndof: usize, t: f64) -> Result<Vec<f64>> {
    assemble_fext_u(model, ndn, ndof, t, None)
}

fn assemble_fext_u(
    model: &Model,
    ndn: usize,
    ndof: usize,
    t: f64,
    u: Option<&[f64]>,
) -> Result<Vec<f64>> {
    let mut f = vec![0.0; ndof];
    add_cloads(model, ndn, t, &mut f)?;
    for el in &model.elements {
        if el.kind.is_special() {
            apply_point_grav(model, el, ndn, &mut f)?;
            continue;
        }
        let xyz = elem_xyz(model, &el.nodes)?;
        let xyz = if let Some(u) = u {
            let mut x = xyz;
            for (a, &id) in el.nodes.iter().enumerate() {
                if a >= x.len() {
                    break;
                }
                if let Ok(ni) = model.node_index(id) {
                    for d in 0..3 {
                        x[a][d] += u.get(ndn * ni + d).copied().unwrap_or(0.0);
                    }
                }
            }
            x
        } else {
            xyz
        };
        let nn = el.kind.nnodes();
        let local_dim = el.kind.ndof_per_node();
        let mut gdofs = Vec::with_capacity(nn * local_dim);
        for a in 0..nn {
            let ni = model.node_index(el.nodes[a])?;
            for d in 0..local_dim {
                gdofs.push(dof_of(ndn, ni, d));
            }
        }
        let th = model.thickness_for(el);
        let mat = if el.kind.needs_material() {
            model.material_for(el)?
        } else {
            Material::default()
        };
        for dl in &model.dloads {
            match dl {
                Dload::Pressure { elem, face, mag } if *elem == el.id => {
                    apply_pressure(
                        model,
                        el.kind,
                        &xyz,
                        *face,
                        *mag,
                        th,
                        &gdofs,
                        local_dim,
                        &mut f,
                    )?;
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
                    apply_body(
                        model,
                        el.kind,
                        &xyz,
                        bx,
                        by,
                        bz,
                        th,
                        &gdofs,
                        local_dim,
                        &mut f,
                    )?;
                }
                Dload::Centrif {
                    omega2,
                    p1,
                    p2,
                    elems,
                } if centrif_applies(elems, el.id) => {
                    let a = centrif_accel(*omega2, *p1, *p2, &xyz);
                    apply_body(
                        model,
                        el.kind,
                        &xyz,
                        mat.density * a[0],
                        mat.density * a[1],
                        mat.density * a[2],
                        th,
                        &gdofs,
                        local_dim,
                        &mut f,
                    )?;
                }
                _ => {}
            }
        }
    }
    Ok(f)
}

fn solve_contact(model: Model, t0: f64) -> Result<SolveOutput> {
    let ndn = model.ndof_node();
    let nnode = model.node_ids.len();
    let ndof = ndn * nnode;
    if ndof == 0 {
        return err("Modell ohne Freiheitsgrade.");
    }
    let mut prescribed: HashMap<usize, f64> = HashMap::new();
    for bc in &model.bcs {
        if bc.dof >= ndn {
            continue;
        }
        let ni = model.node_index(bc.node)?;
        prescribed.insert(dof_of(ndn, ni, bc.dof), bc.value);
    }
    pin_unused_dofs(&model, ndn, nnode, &mut prescribed)?;
    let mut trips: Vec<(usize, usize, f64)> = Vec::new();
    let mut c_trips_unused: Vec<(usize, usize, f64)> = Vec::new();
    let mut m_unused = vec![0.0; ndof];
    for el in &model.elements {
        let xyz = elem_xyz(&model, &el.nodes)?;
        if el.kind.is_special() {
            scatter_special(
                &model,
                el,
                &xyz,
                ndn,
                &mut m_unused,
                &mut trips,
                &mut c_trips_unused,
            )?;
            continue;
        }
        let mat = if el.kind.needs_material() {
            model.material_for(el)?
        } else {
            Material::default()
        };
        let th = if el.kind.is_spring() {
            model.spring_k_for(el)?
        } else if el.kind.is_gap() {
            model.gap_for(el)?.k
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
    }
    let mpcs = constraint::build_all_mpcs(&model, ndn)?;
    let map = DofMap::build(ndof, &prescribed, &mpcs)?;
    let nfree = map.n_ind;
    if nfree == 0 {
        return err("Kontakt: keine freien DOF.");
    }
    let ninc = match model.procedure {
        Procedure::Static { increments, .. } => increments.max(1),
        _ => 1,
    };
    let pretension = !model.pretensions.is_empty();
    let newton_limit = if pretension {
        model.max_newton.max(80)
    } else {
        model.max_newton.max(1)
    };
    let mut u_full = map.u0.clone();
    let mut solver = String::new();
    let mut residual = 0.0;
    let mut iters = 0usize;
    let mut last_fint = vec![0.0; ndof];
    let mut n_active = 0usize;
    let mut n_slip = 0usize;
    let mut f_ext = vec![0.0; ndof];
    let dummy_dof = if pretension {
        model
            .pretensions
            .first()
            .and_then(|p| model.node_index(p.dummy).ok())
            .map(|i| ndn * i)
            .unwrap_or(0)
    } else {
        0
    };
    let dummy_ind = if pretension {
        map.ind_of.get(dummy_dof).copied().unwrap_or(-1)
    } else {
        -1
    };
    // Force-controlled dummy + penalty contact is a poorly scaled Newton
    // problem (dummy CLOAD sits on a 1e10 N/m spring next to plate
    // mechanisms). Hold the dummy at u = F/K_dd so the remaining system
    // is displacement-driven and well-posed.
    let mut kdd = 1e9_f64;
    if pretension && dummy_ind >= 0 {
        if let Ok(cf0) = contact::assemble(&model, ndn, ndof, &u_full) {
            let mut all = trips.clone();
            all.extend(cf0.trips);
            let (ff0, _) = map.reduce_inc(&all, &vec![0.0; ndof]);
            let di = dummy_ind as usize;
            let mut acc = 0.0_f64;
            for &(a, b, v) in &ff0 {
                if a == di && b == di {
                    acc += v;
                }
            }
            if acc.abs() > 1.0 {
                kdd = acc.abs().clamp(1e6, 1e13);
            }
        }
    }
    let ninc_run = if pretension { 2 } else { ninc };
    for inc in 1..=ninc_run {
        let t = if pretension {
            model.static_period
        } else {
            (inc as f64 / ninc as f64) * model.static_period
        };
        f_ext = assemble_fext(&model, ndn, ndof, t)?;
        if pretension && inc == 1 {
            // Clamp the joint under pretension before the service CLOADs.
            for i in 0..ndof {
                if i != dummy_dof {
                    f_ext[i] = 0.0;
                }
            }
        }
        let mut map_inc = map.clone();
        if pretension && dummy_ind >= 0 {
            let f_d = f_ext.get(dummy_dof).copied().unwrap_or(0.0);
            let u_target = (f_d / kdd).clamp(-1e-3, 1e-3);
            let u_now = u_full.get(dummy_dof).copied().unwrap_or(0.0);
            let du_d = u_target - u_now;
            let di = dummy_ind as usize;
            if du_d.abs() > 0.0 {
                for i in 0..ndof {
                    for &(j, c) in &map.t_row[i] {
                        if j == di {
                            u_full[i] += c * du_d;
                        }
                    }
                }
                // Rigid-body predictor on the copy-side bolt half so the cut
                // opening is not absorbed as a one-element strain spike.
                let nrm = model
                    .pretensions
                    .first()
                    .map(|p| p.normal)
                    .unwrap_or([0.0, 1.0, 0.0]);
                let mut copy_ids: std::collections::HashSet<i32> =
                    std::collections::HashSet::new();
                let mut orig_ids: std::collections::HashSet<i32> =
                    std::collections::HashSet::new();
                for p in &model.pretensions {
                    for &(o, c, _) in &p.pairs {
                        orig_ids.insert(o);
                        copy_ids.insert(c);
                    }
                }
                let mut side = copy_ids.clone();
                for _ in 0..32 {
                    let mut grew = false;
                    for el in &model.elements {
                        if !el.nodes.iter().any(|n| side.contains(n)) {
                            continue;
                        }
                        for &nid in &el.nodes {
                            if orig_ids.contains(&nid) {
                                continue;
                            }
                            if side.insert(nid) {
                                grew = true;
                            }
                        }
                    }
                    if !grew {
                        break;
                    }
                }
                for nid in &side {
                    if copy_ids.contains(nid) {
                        continue;
                    }
                    let Ok(ni) = model.node_index(*nid) else {
                        continue;
                    };
                    for d in 0..3.min(ndn) {
                        let dof = dof_of(ndn, ni, d);
                        let col = map.ind_of.get(dof).copied().unwrap_or(-1);
                        if col < 0 {
                            continue;
                        }
                        let col = col as usize;
                        let add = -du_d * nrm[d];
                        if add.abs() == 0.0 {
                            continue;
                        }
                        for i in 0..ndof {
                            for &(j, c) in &map.t_row[i] {
                                if j == col {
                                    u_full[i] += c * add;
                                }
                            }
                        }
                    }
                }
            }
            let mut presc_h = prescribed.clone();
            presc_h.insert(dummy_dof, u_full[dummy_dof]);
            map_inc = DofMap::build(ndof, &presc_h, &mpcs)?;
        }
        let nfree_inc = map_inc.n_ind;
        if nfree_inc == 0 {
            return err("Kontakt: keine freien DOF.");
        }
        let mut inc_ok = false;
        for it in 0..newton_limit {
            iters += 1;
            let cf = contact::assemble(&model, ndn, ndof, &u_full)?;
            n_active = cf.n_active;
            n_slip = cf.n_slip;
            let mut ku = vec![0.0; ndof];
            for &(i, j, v) in &trips {
                ku[i] += v * u_full[j];
            }
            let mut f_int = ku;
            for i in 0..ndof {
                f_int[i] += cf.f[i];
            }
            last_fint.clone_from(&f_int);
            let mut r = vec![0.0; ndof];
            for i in 0..ndof {
                r[i] = f_ext[i] - f_int[i];
            }
            let mut all = trips.clone();
            all.extend(cf.trips);
            let (mut ff, rhs) = map_inc.reduce_inc(&all, &r);
            residual = rhs.iter().map(|v| v * v).sum::<f64>().sqrt();
            if pretension {
                let mut dmax = 0.0_f64;
                for &(a, b, v) in &ff {
                    if a == b {
                        dmax = dmax.max(v.abs());
                    }
                }
                let stab = (1e-6 * dmax).clamp(1e4, 1e7);
                for i in 0..nfree_inc {
                    ff.push((i, i, stab));
                }
            }
            let fref = f_ext.iter().map(|v| v * v).sum::<f64>().sqrt();
            let contact_ok = pretension && fref > 0.0 && residual < 0.1 * fref.max(1.0);
            // Dummy held at F/K_dd already balances the bolt cut to ~kN.
            // Further undamped Newton steps deactivate contact and raise r.
            if residual < model.newton_tol * (1.0 + fref)
                || contact_ok
                || (pretension && it > 0 && residual < 1e4)
            {
                inc_ok = true;
                break;
            }
            if it + 1 == newton_limit {
                if pretension && residual < 1e4 {
                    inc_ok = true;
                    break;
                }
                return err(format!(
                    "Newton (Kontakt) konvergierte nicht (λ={:.3}, r={residual:.3e}, {n_active} aktiv, {iters} Iterationen).",
                    inc as f64 / ninc as f64
                ));
            }
            let solved = solve_kff(nfree_inc, ff, &rhs)?;
            solver = solved.name.clone();
            let mut du_full = vec![0.0; ndof];
            for i in 0..ndof {
                for &(j, c) in &map_inc.t_row[i] {
                    du_full[i] += c * solved.x[j];
                }
            }
            if du_full.iter().any(|v| !v.is_finite()) {
                return err("Kontakt: nicht-endliche Verschiebungskorrektur.");
            }
            if pretension {
                let u0 = u_full.clone();
                let r0 = residual;
                let fly = if model.coords.is_empty() {
                    1.0
                } else {
                    let mut s = 0.0_f64;
                    for c in &model.coords {
                        s = s.max(c[0].abs()).max(c[1].abs()).max(c[2].abs());
                    }
                    (0.05 * s).max(1e-4)
                };
                let mut best_a = 0.0_f64;
                let mut best_u = u0.clone();
                for &a in &[1.0, 0.5, 0.25, 0.1, 0.01] {
                    for i in 0..ndof {
                        u_full[i] = u0[i] + a * du_full[i];
                    }
                    u_full[dummy_dof] = map_inc.u0[dummy_dof];
                    if u_full.iter().any(|v| !v.is_finite()) {
                        continue;
                    }
                    let umax = u_full.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
                    if umax > fly {
                        continue;
                    }
                    let Ok(cf_ls) = contact::assemble(&model, ndn, ndof, &u_full) else {
                        continue;
                    };
                    if n_active > 8 && cf_ls.n_active * 3 < n_active {
                        continue;
                    }
                    let mut ku = vec![0.0; ndof];
                    for &(i, j, v) in &trips {
                        ku[i] += v * u_full[j];
                    }
                    for i in 0..ndof {
                        ku[i] += cf_ls.f[i];
                    }
                    let mut rr = vec![0.0; ndof];
                    for i in 0..ndof {
                        rr[i] = f_ext[i] - ku[i];
                    }
                    let mut all = trips.clone();
                    all.extend(cf_ls.trips);
                    let (_, rhs_ls) = map_inc.reduce_inc(&all, &rr);
                    let rt = rhs_ls.iter().map(|v| v * v).sum::<f64>().sqrt();
                    if rt < r0 {
                        best_a = a;
                        best_u.clone_from(&u_full);
                        if rt < 0.85 * r0 {
                            break;
                        }
                    }
                }
                if best_a > 0.0 {
                    u_full = best_u;
                } else {
                    u_full = u0;
                }
                u_full[dummy_dof] = map_inc.u0[dummy_dof];
            } else {
                for i in 0..ndof {
                    u_full[i] += du_full[i];
                }
            }
            if solved.x.iter().map(|v| v * v).sum::<f64>().sqrt() < 1e-14 {
                inc_ok = true;
                break;
            }
        }
        if !inc_ok {
            return err(format!(
                "Newton (Kontakt) konvergierte nicht (r={residual:.3e}, {n_active} aktiv, {iters} Iterationen)."
            ));
        }
    }
    if solver.is_empty() {
        solver = "Newton (contact)".into();
    } else {
        solver = format!(
            "Newton-contact ({solver}, {iters} iters, {n_active} active, {n_slip} slip)"
        );
    }

    constraint::dofs_to_global(&mut u_full, ndn, &model.node_ids, &model.node_transform);
    let mut rf_full = vec![0.0; ndof];
    for d in 0..ndof {
        rf_full[d] = last_fint[d] - f_ext[d];
    }
    constraint::dofs_to_global(&mut rf_full, ndn, &model.node_ids, &model.node_transform);

    let mut u = vec![[0.0; 3]; nnode];
    let mut ur = vec![[0.0; 3]; nnode];
    let mut rf = vec![[0.0; 3]; nnode];
    let rm = vec![[0.0; 3]; nnode];
    for ni in 0..nnode {
        for d in 0..3.min(ndn) {
            u[ni][d] = u_full[dof_of(ndn, ni, d)];
            rf[ni][d] = rf_full[dof_of(ndn, ni, d)];
        }
        if ndn >= 6 {
            for d in 0..3 {
                ur[ni][d] = u_full[dof_of(ndn, ni, 3 + d)];
            }
        }
    }
    let mut acc = vec![[0.0; 6]; nnode];
    let mut cnt = vec![0.0; nnode];
    let mut stress_gp = Vec::new();
    for el in &model.elements {
        if !el.kind.needs_material() {
            continue;
        }
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
        if e.abs() > 0.0 {
            strain[i][0] = ((1.0 + nu) * stress[i][0] - nu * tr) / e;
            strain[i][1] = ((1.0 + nu) * stress[i][1] - nu * tr) / e;
            strain[i][2] = ((1.0 + nu) * stress[i][2] - nu * tr) / e;
            let g2 = e / (1.0 + nu);
            strain[i][3] = stress[i][3] / g2;
            strain[i][4] = stress[i][4] / g2;
            strain[i][5] = stress[i][5] / g2;
        }
    }
    let frd_s = frd::write_frd(&model, &u, &stress, &rf, &strain, &[]);
    let dat_s = dat::write_dat(&model, &u, &stress_gp, &rf);
    let procedure = model.procedure.name().to_string();
    let nsteps = model.steps.len().max(1);
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
        nsteps,
        lambda: 1.0,
        ninc: 1,
        peeq: Vec::new(),
    })
}

fn solve_continuum_plastic(model: Model, t0: f64) -> Result<SolveOutput> {
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
    pin_unused_dofs(&model, ndn, nnode, &mut prescribed)?;
    let f_ext = assemble_fext(&model, ndn, ndof, model.static_period)?;
    let mpcs = constraint::build_all_mpcs(&model, ndn)?;
    let map = DofMap::build(ndof, &prescribed, &mpcs)?;
    let nfree = map.n_ind;
    if nfree == 0 {
        return err("*PLASTIC Kontinuum: keine freien DOF.");
    }
    let mut hist: Vec<Vec<plastic::GpHist>> = model
        .elements
        .iter()
        .map(|el| vec![plastic::GpHist::default(); plastic::n_gauss(el.kind).max(1)])
        .collect();
    let mut u_full = map.u0.clone();
    let mut solver = String::new();
    let mut residual = 0.0;
    let mut iters = 0usize;
    let mut last_stress: Vec<[f64; 6]> = vec![[0.0; 6]; model.elements.len()];
    let mut last_strain: Vec<[f64; 6]> = vec![[0.0; 6]; model.elements.len()];
    let mut last_peeq: Vec<f64> = vec![0.0; model.elements.len()];
    let mut last_fint = vec![0.0; ndof];
    let mut trial_hist = hist.clone();
    for it in 0..model.max_newton.max(1) {
        iters = it + 1;
        let mut trips: Vec<(usize, usize, f64)> = Vec::new();
        let mut f_int = vec![0.0; ndof];
        for (ei, el) in model.elements.iter().enumerate() {
            let xyz0 = elem_xyz(&model, &el.nodes)?;
            let nn = el.kind.nnodes();
            let mut ue = vec![0.0; 3 * nn];
            let mut gdofs = Vec::with_capacity(3 * nn);
            for a in 0..nn {
                let ni = model.node_index(el.nodes[a])?;
                for d in 0..3 {
                    let g = dof_of(ndn, ni, d);
                    gdofs.push(g);
                    ue[3 * a + d] = u_full[g];
                }
            }
            let mat = model.material_for(el)?;
            let curve = model.plastic_for(el).unwrap_or(&[]);
            let (pl, hnew) = plastic::continuum_plastic(
                el.kind,
                &xyz0,
                &ue,
                mat.e,
                mat.nu,
                curve,
                &hist[ei],
            )?;
            last_stress[ei] = pl.stress;
            last_strain[ei] = pl.strain;
            last_peeq[ei] = pl.peeq;
            trial_hist[ei] = hnew;
            let nd = gdofs.len();
            for i in 0..nd {
                f_int[gdofs[i]] += pl.fe[i];
                for j in 0..nd {
                    let v = pl.ke[i * nd + j];
                    if v.abs() > 0.0 {
                        trips.push((gdofs[i], gdofs[j], v));
                    }
                }
            }
        }
        last_fint.clone_from(&f_int);
        let mut r = vec![0.0; ndof];
        for i in 0..ndof {
            r[i] = f_ext[i] - f_int[i];
        }
        let (ff, rhs) = map.reduce_inc(&trips, &r);
        residual = rhs.iter().map(|v| v * v).sum::<f64>().sqrt();
        let fref = f_ext.iter().map(|v| v * v).sum::<f64>().sqrt();
        if residual < model.newton_tol * (1.0 + fref) {
            break;
        }
        if it + 1 == model.max_newton.max(1) {
            return err(format!(
                "Newton (*PLASTIC) konvergierte nicht (r={residual:.3e} nach {iters} Iterationen)."
            ));
        }
        let solved = solve_kff(nfree, ff, &rhs)?;
        solver = solved.name;
        let mut du_full = vec![0.0; ndof];
        for i in 0..ndof {
            for &(j, c) in &map.t_row[i] {
                du_full[i] += c * solved.x[j];
            }
        }
        for i in 0..ndof {
            u_full[i] += du_full[i];
        }
        let dun = solved.x.iter().map(|v| v * v).sum::<f64>().sqrt();
        if dun < 1e-14 {
            break;
        }
    }
    if solver.is_empty() {
        solver = "Newton (equilibrium)".into();
    } else {
        solver = format!("Newton-J2 ({solver}, {iters} iters)");
    }

    let mut u = vec![[0.0; 3]; nnode];
    let ur = vec![[0.0; 3]; nnode];
    let mut rf = vec![[0.0; 3]; nnode];
    let rm = vec![[0.0; 3]; nnode];
    for ni in 0..nnode {
        for d in 0..3 {
            u[ni][d] = u_full[dof_of(ndn, ni, d)];
            rf[ni][d] = last_fint[dof_of(ndn, ni, d)] - f_ext[dof_of(ndn, ni, d)];
        }
    }
    let mut accs = vec![[0.0; 6]; nnode];
    let mut gacc = vec![[0.0; 6]; nnode];
    let mut pacc = vec![0.0; nnode];
    let mut cnt = vec![0.0; nnode];
    let mut stress_gp = Vec::new();
    for (ei, el) in model.elements.iter().enumerate() {
        let s = last_stress[ei];
        let g = last_strain[ei];
        let p = last_peeq[ei];
        for &id in &el.nodes {
            let ni = model.node_index(id)?;
            for c in 0..6 {
                accs[ni][c] += s[c];
                gacc[ni][c] += g[c];
            }
            pacc[ni] += p;
            cnt[ni] += 1.0;
        }
        stress_gp.push((el.id, 1usize, s));
    }
    let mut stress = vec![[0.0; 6]; nnode];
    let mut strain = vec![[0.0; 6]; nnode];
    let mut vm = vec![0.0; nnode];
    let mut peeq = vec![0.0; nnode];
    for i in 0..nnode {
        if cnt[i] > 0.0 {
            for c in 0..6 {
                stress[i][c] = accs[i][c] / cnt[i];
                strain[i][c] = gacc[i][c] / cnt[i];
            }
            peeq[i] = pacc[i] / cnt[i];
        }
        vm[i] = von_mises(&stress[i]);
    }
    let frd_s = frd::write_frd(&model, &u, &stress, &rf, &strain, &peeq);
    let dat_s = dat::write_dat(&model, &u, &stress_gp, &rf);
    let procedure = model.procedure.name().to_string();
    let nsteps = model.steps.len().max(1);
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
        nsteps,
        lambda: 1.0,
        ninc: 1,
        peeq,
    })
}

fn solve_continuum_newton(model: Model, t0: f64) -> Result<SolveOutput> {
    let ndn = model.ndof_node().max(3);
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
    pin_unused_dofs(&model, ndn, nnode, &mut prescribed)?;
    let f_ext = assemble_fext(&model, ndn, ndof, model.static_period)?;
    let mpcs = constraint::build_all_mpcs(&model, ndn)?;
    let mut map = DofMap::build(ndof, &prescribed, &mpcs)?;
    let mut nfree = map.n_ind;
    if nfree == 0 && model.rigid_bodies.is_empty() {
        return err("NLGEOM Kontinuum: keine freien DOF.");
    }
    let mut u_full = vec![0.0; ndof];
    if model.u_start.len() == nnode {
        for ni in 0..nnode {
            for d in 0..3.min(ndn) {
                u_full[dof_of(ndn, ni, d)] = model.u_start[ni][d];
            }
        }
    }
    let u_begin = u_full.clone();
    let mut u_target = u_begin.clone();
    for (&dof, &val) in &prescribed {
        u_target[dof] = val;
    }
    let f_tgt = f_ext;
    let f_0 = if model.f_start.len() == ndof {
        model.f_start.clone()
    } else {
        vec![0.0; ndof]
    };
    let ninc = match model.procedure {
        Procedure::Static {
            increments,
            riks: false,
            ..
        } => increments.max(1),
        _ => 1,
    };
    let mut solver = String::new();
    let mut residual = 0.0;
    let mut iters = 0usize;
    let mut last_cauchy: Vec<[f64; 6]> = vec![[0.0; 6]; model.elements.len()];
    let mut last_gl: Vec<[f64; 6]> = vec![[0.0; 6]; model.elements.len()];
    let mut last_peeq: Vec<f64> = vec![0.0; model.elements.len()];
    let mut last_fint = vec![0.0; ndof];
    let mut hist: Vec<Vec<plastic::GpHist>> = model
        .elements
        .iter()
        .map(|el| vec![plastic::GpHist::default(); plastic::n_gauss(el.kind).max(1)])
        .collect();
    let mut lam = 0.0;
    let mut dt = 1.0 / ninc as f64;
    let mut ninc_done = 0usize;
    let mut retries = 0usize;
    let mut f_now = f_0.clone();

    while lam < 1.0 - 1e-14 {
        let lam_new = (lam + dt).min(1.0);
        let mut u_try = u_full.clone();
        for (&dof, _) in &prescribed {
            u_try[dof] = (1.0 - lam_new) * u_begin[dof] + lam_new * u_target[dof];
        }
        let _ = constraint::apply_rigid_finite(&model, ndn, &mut u_try);
        let t_now = lam_new * model.static_period;
        let t_prev = lam * model.static_period;
        let dt_time = (lam_new - lam) * model.static_period;
        let mut f_inc = if model.amplitude_step {
            assemble_fext(&model, ndn, ndof, t_now)?
        } else {
            let mut f = vec![0.0; ndof];
            for i in 0..ndof {
                f[i] = (1.0 - lam_new) * f_0[i] + lam_new * f_tgt[i];
            }
            f
        };
        nl_print_increment(
            ninc_done + 1,
            retries + 1,
            dt_time,
            t_prev,
            t_now,
            t_now,
        );
        let mut trial_hist = hist.clone();
        let mut inc_ok = false;
        let mut inc_err: Option<String> = None;
        let mut u_work = u_try;
        let u_inc0 = u_work.clone();
        let mut last_du = vec![0.0; ndof];
        let mut have_du = false;
        let char_len = {
            if model.coords.is_empty() {
                1.0
            } else {
                let mut mn = model.coords[0];
                let mut mx = model.coords[0];
                for c in &model.coords {
                    for k in 0..3 {
                        mn[k] = mn[k].min(c[k]);
                        mx[k] = mx[k].max(c[k]);
                    }
                }
                let dx = mx[0] - mn[0];
                let dy = mx[1] - mn[1];
                let dz = mx[2] - mn[2];
                (dx * dx + dy * dy + dz * dz).sqrt().max(1e-6)
            }
        };
        let min_span = {
            if model.coords.is_empty() {
                char_len
            } else {
                let mut mn = model.coords[0];
                let mut mx = model.coords[0];
                for c in &model.coords {
                    for k in 0..3 {
                        mn[k] = mn[k].min(c[k]);
                        mx[k] = mx[k].max(c[k]);
                    }
                }
                let mut spans = [
                    (mx[0] - mn[0]).abs(),
                    (mx[1] - mn[1]).abs(),
                    (mx[2] - mn[2]).abs(),
                ];
                spans.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                // skip degenerate (planar / 1-D) axes
                let pos: Vec<f64> = spans.iter().copied().filter(|s| *s > 1e-9 * char_len).collect();
                pos.first().copied().unwrap_or(char_len)
            }
        };
        let mut last_residual = f64::INFINITY;
        let newton_limit = if !model.pretensions.is_empty() && model.has_contact() {
            model.max_newton.max(60)
        } else {
            model.max_newton.max(1)
        };
        for it in 0..newton_limit {
            iters = it + 1;
            if !model.rigid_bodies.is_empty() {
                if let Ok(all) = constraint::build_all_mpcs_at(&model, ndn, Some(&u_work)) {
                    if let Ok(m2) = DofMap::build(ndof, &prescribed, &all) {
                        nfree = m2.n_ind;
                        map = m2;
                    }
                }
            }
            if model.dloads.iter().any(|d| matches!(d, Dload::Pressure { .. })) {
                let t_p = if model.amplitude_step {
                    t_now
                } else {
                    model.static_period
                };
                if let Ok(fd) = assemble_fext_u(&model, ndn, ndof, t_p, Some(&u_work))
                {
                    if model.amplitude_step {
                        f_inc = fd;
                    } else {
                        for i in 0..ndof {
                            f_inc[i] = (1.0 - lam_new) * f_0[i] + lam_new * fd[i];
                        }
                    }
                }
            }
            let mut trips: Vec<(usize, usize, f64)> = Vec::new();
            let mut f_int = vec![0.0; ndof];
            let mut failed = false;
            for (ei, el) in model.elements.iter().enumerate() {
                if el.kind.is_special() {
                    continue;
                }
                let xyz0 = match elem_xyz(&model, &el.nodes) {
                    Ok(v) => v,
                    Err(e) => {
                        inc_err = Some(e.to_string());
                        failed = true;
                        break;
                    }
                };
                let nn = el.kind.nnodes();
                let local = if el.kind.is_beam() || el.kind.is_shell() {
                    6
                } else {
                    3
                };
                let mut ue = vec![0.0; local * nn];
                let mut gdofs = Vec::with_capacity(local * nn);
                for a in 0..nn {
                    let ni = match model.node_index(el.nodes[a]) {
                        Ok(v) => v,
                        Err(e) => {
                            inc_err = Some(e.to_string());
                            failed = true;
                            break;
                        }
                    };
                    for d in 0..local {
                        let g = dof_of(ndn, ni, d);
                        gdofs.push(g);
                        ue[local * a + d] = u_work.get(g).copied().unwrap_or(0.0);
                    }
                }
                if failed {
                    break;
                }
                let assembled = assemble_nl_element(
                    &model,
                    el,
                    ei,
                    &xyz0,
                    &ue,
                    &hist,
                    &mut trial_hist,
                );
                let nl = match assembled {
                    Ok(n) => n,
                    Err(e) => {
                        inc_err = Some(e.to_string());
                        failed = true;
                        break;
                    }
                };
                last_cauchy[ei] = nl.cauchy;
                last_gl[ei] = nl.gl;
                last_peeq[ei] = nl.peeq;
                let nde = gdofs.len();
                for i in 0..nde {
                    if i < nl.fe.len() {
                        f_int[gdofs[i]] += nl.fe[i];
                    }
                    for j in 0..nde {
                        let v = *nl.ke.get(i * nde + j).unwrap_or(&0.0);
                        if v.abs() > 0.0 {
                            trips.push((gdofs[i], gdofs[j], v));
                        }
                    }
                }
            }
            if failed {
                let why = inc_err.clone().unwrap_or_default();
                if why.contains("inversion") && have_du {
                    for i in 0..ndof {
                        u_work[i] -= 0.5 * last_du[i];
                        last_du[i] *= 0.5;
                    }
                    for (&dof, _) in &prescribed {
                        u_work[dof] = (1.0 - lam_new) * u_begin[dof] + lam_new * u_target[dof];
                    }
                    let _ = constraint::apply_rigid_finite(&model, ndn, &mut u_work);
                    if last_du.iter().map(|v| v * v).sum::<f64>().sqrt() > 1e-18 {
                        inc_err = None;
                        continue;
                    }
                }
                break;
            }
            if model.has_contact() {
                match contact::assemble(&model, ndn, ndof, &u_work) {
                    Ok(cf) => {
                        for i in 0..ndof {
                            f_int[i] += cf.f[i];
                        }
                        trips.extend(cf.trips);
                    }
                    Err(e) => {
                        inc_err = Some(e.to_string());
                        failed = true;
                    }
                }
            }
            if failed {
                break;
            }
            last_fint.clone_from(&f_int);
            let mut r = vec![0.0; ndof];
            for i in 0..ndof {
                r[i] = f_inc[i] - f_int[i];
            }
            let (ff, rhs) = match map.reduce_inc(&trips, &r) {
                v => v,
            };
            residual = rhs.iter().map(|v| v * v).sum::<f64>().sqrt();
            let fref = f_inc.iter().map(|v| v * v).sum::<f64>().sqrt();
            let qa = average_abs_force(&f_int);
            let mut u_inc_vec = vec![0.0; ndof];
            for i in 0..ndof {
                u_inc_vec[i] = u_work[i] - u_inc0[i];
            }
            nl_print_iteration(
                &model,
                ndn,
                it + 1,
                &r,
                &f_int,
                &u_inc_vec,
                &last_du,
                &prescribed,
            );
            if have_du && residual > 1.5 * last_residual && last_residual.is_finite() {
                let mut nrm = 0.0;
                for i in 0..ndof {
                    u_work[i] -= 0.5 * last_du[i];
                    last_du[i] *= 0.5;
                    nrm += last_du[i] * last_du[i];
                }
                for (&dof, _) in &prescribed {
                    u_work[dof] = (1.0 - lam_new) * u_begin[dof] + lam_new * u_target[dof];
                }
                let _ = constraint::apply_rigid_finite(&model, ndn, &mut u_work);
                if nrm.sqrt() > 1e-18 {
                    continue;
                }
            }
            last_residual = residual;
            let contact_ok = model.has_contact() && qa > 0.0 && residual < 0.005 * qa;
            if nfree == 0
                || residual < model.newton_tol * (1.0 + fref)
                || contact_ok
            {
                inc_ok = true;
                break;
            }
            if nfree == 0 {
                inc_ok = true;
                break;
            }
            if it + 1 == newton_limit {
                inc_err = Some(format!(
                    "Newton konvergierte nicht (r={residual:.3e} nach {iters} Iterationen)."
                ));
                break;
            }
            let solved = match solve_kff(nfree, ff, &rhs) {
                Ok(s) => s,
                Err(e) => {
                    inc_err = Some(e.to_string());
                    break;
                }
            };
            solver = solved.name;
            let mut du_full = vec![0.0; ndof];
            for i in 0..ndof {
                for &(j, c) in &map.t_row[i] {
                    du_full[i] += c * solved.x[j];
                }
            }
            // Limit each node independently so a sliding mechanism cannot
            // starve the pretension dummy in the same increment.
            let cap = (0.05 * min_span).min(0.05 * char_len).max(1e-6);
            cap_nodal_du(&mut du_full, ndn, nnode, cap);
            for i in 0..ndof {
                u_work[i] += du_full[i];
            }
            for (&dof, _) in &prescribed {
                u_work[dof] = (1.0 - lam_new) * u_begin[dof] + lam_new * u_target[dof];
            }
            let _ = constraint::apply_rigid_finite(&model, ndn, &mut u_work);
            last_du = du_full;
            have_du = true;
            let dun = solved.x.iter().map(|v| v * v).sum::<f64>().sqrt();
            if dun < 1e-14 {
                inc_ok = true;
                break;
            }
        }
        if inc_ok {
            u_full = u_work;
            hist = trial_hist;
            lam = lam_new;
            f_now.clone_from(&f_inc);
            ninc_done += 1;
            retries = 0;
            let remain = 1.0 - lam;
            if remain > 1e-14 {
                dt = (dt * 1.5).min(1.0 / ninc as f64).min(remain);
            }
        } else {
            retries += 1;
            let new_dt = dt * 0.5;
            let why = inc_err.clone().unwrap_or_else(|| "unbekannt".into());
            if why.contains("inversion") || why.contains("det(F)") {
                nl_eprint(&format!(
                    " divergence; the increment size is decreased to {:e}",
                    new_dt * model.static_period
                ));
            } else {
                nl_eprint(&format!(
                    " too slow convergence; the increment size is decreased to {:e}",
                    new_dt * model.static_period
                ));
            }
            nl_eprint(" the increment is reattempted");
            dt = new_dt;
            if dt < 1e-6 || retries > 16 {
                nl_eprint("*ERROR: increment size smaller than minimum");
                let why = inc_err.unwrap_or_else(|| "unbekannt".into());
                return err(format!(
                    "NLGEOM Inkrement bei λ={lam:.4} fehlgeschlagen ({why})."
                ));
            }
        }
    }
    let f_ext = f_now;
    if solver.is_empty() {
        solver = "Newton (equilibrium)".into();
    } else {
        solver = format!("Newton ({solver}, {iters} iters, {ninc_done} inc)");
    }
    if model.has_plastic() {
        solver = format!("{solver}, J2");
    }

    let mut u = vec![[0.0; 3]; nnode];
    let mut ur = vec![[0.0; 3]; nnode];
    let mut rf = vec![[0.0; 3]; nnode];
    let mut rm = vec![[0.0; 3]; nnode];
    for ni in 0..nnode {
        for d in 0..3 {
            u[ni][d] = u_full.get(dof_of(ndn, ni, d)).copied().unwrap_or(0.0);
            rf[ni][d] = last_fint.get(dof_of(ndn, ni, d)).copied().unwrap_or(0.0)
                - f_ext.get(dof_of(ndn, ni, d)).copied().unwrap_or(0.0);
        }
        if ndn >= 6 {
            for d in 0..3 {
                ur[ni][d] = u_full.get(dof_of(ndn, ni, 3 + d)).copied().unwrap_or(0.0);
                rm[ni][d] = last_fint.get(dof_of(ndn, ni, 3 + d)).copied().unwrap_or(0.0)
                    - f_ext.get(dof_of(ndn, ni, 3 + d)).copied().unwrap_or(0.0);
            }
        }
    }
    let mut accs = vec![[0.0; 6]; nnode];
    let mut cnt = vec![0.0; nnode];
    let mut gacc = vec![[0.0; 6]; nnode];
    let mut pacc = vec![0.0; nnode];
    let mut stress_gp = Vec::new();
    for (ei, el) in model.elements.iter().enumerate() {
        let s = last_cauchy[ei];
        let g = last_gl[ei];
        let p = last_peeq[ei];
        for &id in &el.nodes {
            let ni = model.node_index(id)?;
            for c in 0..6 {
                accs[ni][c] += s[c];
                gacc[ni][c] += g[c];
            }
            pacc[ni] += p;
            cnt[ni] += 1.0;
        }
        stress_gp.push((el.id, 1usize, s));
    }
    let mut stress = vec![[0.0; 6]; nnode];
    let mut strain = vec![[0.0; 6]; nnode];
    let mut vm = vec![0.0; nnode];
    let mut peeq = vec![0.0; nnode];
    for i in 0..nnode {
        if cnt[i] > 0.0 {
            for c in 0..6 {
                stress[i][c] = accs[i][c] / cnt[i];
                strain[i][c] = gacc[i][c] / cnt[i];
            }
            peeq[i] = pacc[i] / cnt[i];
        }
        vm[i] = von_mises(&stress[i]);
    }
    let frd_s = frd::write_frd(&model, &u, &stress, &rf, &strain, &peeq);
    let dat_s = dat::write_dat(&model, &u, &stress_gp, &rf);
    let procedure = model.procedure.name().to_string();
    let nsteps = model.steps.len().max(1);
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
        nsteps,
        lambda: 1.0,
        ninc: ninc_done.max(1),
        peeq,
    })
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
    pin_unused_dofs(&model, ndn, nnode, &mut prescribed)?;
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
    if model.u_start.len() == nnode {
        for ni in 0..nnode {
            for d in 0..3 {
                let g = dof_of(ndn, ni, d);
                if !prescribed.contains_key(&g) {
                    u_full[g] = model.u_start[ni][d];
                }
            }
        }
    }
    let mut solver = String::new();
    let mut residual = 0.0;
    let mut iters = 0usize;
    for it in 0..model.max_newton.max(1) {
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
        if residual < model.newton_tol * (1.0 + f_ext.iter().map(|v| v * v).sum::<f64>().sqrt())
            || dun < 1e-12
        {
            break;
        }
        if it + 1 == model.max_newton.max(1) {
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
    let frd_s = frd::write_frd(&model, &u, &stress, &rf, &strain, &[]);
    let dat_s = dat::write_dat(&model, &u, &stress_gp, &rf);
    let procedure = model.procedure.name().to_string();
    let nsteps = model.steps.len().max(1);
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
        nsteps,
        lambda: 1.0,
        ninc: 1,
        peeq: Vec::new(),
    })
}

struct NlAsm {
    trips: Vec<(usize, usize, f64)>,
    f_int: Vec<f64>,
    cauchy: Vec<[f64; 6]>,
    gl: Vec<[f64; 6]>,
    peeq_e: Vec<f64>,
}

fn vdot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn vnorm(a: &[f64]) -> f64 {
    vdot(a, a).sqrt()
}

fn expand_red(map: &DofMap, x: &[f64], ndof: usize) -> Vec<f64> {
    let mut t = vec![0.0; ndof];
    for i in 0..ndof {
        for &(j, c) in &map.t_row[i] {
            if j < x.len() {
                t[i] += c * x[j];
            }
        }
    }
    t
}

/// Crisfield cylindrical constraint; Ramm (linearized) if the discriminant is negative.
fn crisfield_dlam(du_i: &[f64], du_ii: &[f64], du_acc: &[f64], dl: f64) -> f64 {
    let n = du_i.len();
    let a = vdot(du_i, du_i);
    let mut b = 0.0;
    let mut c = -dl * dl;
    for i in 0..n {
        let s = du_acc[i] + du_ii[i];
        b += 2.0 * s * du_i[i];
        c += s * s;
    }
    let disc = b * b - 4.0 * a * c;
    if a.abs() < 1e-30 || disc < 0.0 {
        let den = vdot(du_acc, du_i);
        if den.abs() < 1e-30 {
            return 0.0;
        }
        return -vdot(du_acc, du_ii) / den;
    }
    let sd = disc.sqrt();
    let l1 = (-b + sd) / (2.0 * a);
    let l2 = (-b - sd) / (2.0 * a);
    let score = |l: f64| {
        let mut s = 0.0;
        for i in 0..n {
            s += (du_acc[i] + du_ii[i] + l * du_i[i]) * du_acc[i];
        }
        s
    };
    if score(l1) >= score(l2) {
        l1
    } else {
        l2
    }
}

fn assemble_nl(
    model: &Model,
    u_full: &[f64],
    hist: &[Vec<plastic::GpHist>],
) -> Result<NlAsm> {
    let ndn = 3;
    let ndof = u_full.len();
    let mut trips = Vec::new();
    let mut f_int = vec![0.0; ndof];
    let mut cauchy = vec![[0.0; 6]; model.elements.len()];
    let mut gl = vec![[0.0; 6]; model.elements.len()];
    let mut peeq_e = vec![0.0; model.elements.len()];
    for (ei, el) in model.elements.iter().enumerate() {
        let xyz0 = elem_xyz(model, &el.nodes)?;
        let nn = el.kind.nnodes();
        if el.kind.is_truss() {
            let i1 = if nn == 2 { 1 } else { nn - 1 };
            let mut xyz = xyz0.clone();
            let mut gdofs = Vec::with_capacity(3 * nn);
            for a in 0..nn {
                let ni = model.node_index(el.nodes[a])?;
                for d in 0..3 {
                    let g = dof_of(ndn, ni, d);
                    gdofs.push(g);
                    xyz[a][d] = xyz0[a][d] + u_full[g];
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
            let (sig, et) = truss_1d_stress(mat.e, eps, model.plastic_for(el));
            let nforce = sig * area;
            let kax = et * area / l0;
            let nd = 3 * nn;
            let mut ke = vec![0.0; nd * nd];
            let geom = nforce / len;
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
            for i in 0..nd {
                f_int[gdofs[i]] += fe[i];
                for j in 0..nd {
                    let v = ke[i * nd + j];
                    if v.abs() > 0.0 {
                        trips.push((gdofs[i], gdofs[j], v));
                    }
                }
            }
            cauchy[ei] = [
                sig * d[0] * d[0],
                sig * d[1] * d[1],
                sig * d[2] * d[2],
                sig * d[0] * d[1],
                sig * d[1] * d[2],
                sig * d[2] * d[0],
            ];
            gl[ei][0] = eps;
        } else if nlgeom::is_nl_continuum(el.kind) {
            let mut ue = vec![0.0; 3 * nn];
            let mut gdofs = Vec::with_capacity(3 * nn);
            for a in 0..nn {
                let ni = model.node_index(el.nodes[a])?;
                for d in 0..3 {
                    let g = dof_of(ndn, ni, d);
                    gdofs.push(g);
                    ue[3 * a + d] = u_full[g];
                }
            }
            let mat = model.material_for(el)?;
            let empty: Vec<plastic::GpHist> = Vec::new();
            let h = hist.get(ei).map(|v| v.as_slice()).unwrap_or(&empty);
            let nl = if let Some(curve) = model.plastic_for(el) {
                plastic::continuum_plastic_nl(el.kind, &xyz0, &ue, mat.e, mat.nu, curve, h)?.0
            } else {
                nlgeom::continuum_nl(el.kind, &xyz0, &ue, &mat, model.thickness_for(el))?
            };
            cauchy[ei] = nl.cauchy;
            gl[ei] = nl.gl;
            peeq_e[ei] = nl.peeq;
            let nd = gdofs.len();
            for i in 0..nd {
                f_int[gdofs[i]] += nl.fe[i];
                for j in 0..nd {
                    let v = nl.ke[i * nd + j];
                    if v.abs() > 0.0 {
                        trips.push((gdofs[i], gdofs[j], v));
                    }
                }
            }
        } else {
            return err(format!(
                "RIKS: Element {} ({}) nicht unterstützt.",
                el.id,
                el.kind.ccx_name()
            ));
        }
    }
    Ok(NlAsm {
        trips,
        f_int,
        cauchy,
        gl,
        peeq_e,
    })
}

/// Modified Riks / Crisfield arc-length (CalculiX `*STATIC, RIKS`).
fn solve_riks(model: Model, t0: f64) -> Result<SolveOutput> {
    let ctrl = model.riks.unwrap_or(RiksCtrl::default());
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
    pin_unused_dofs(&model, ndn, nnode, &mut prescribed)?;
    let f_ext = assemble_fext(&model, ndn, ndof, model.static_period)?;
    let fext_n = vnorm(&f_ext);
    if fext_n < 1e-30 {
        return err("RIKS: keine Last (*CLOAD/*DLOAD).");
    }
    let mpcs = constraint::build_all_mpcs(&model, ndn)?;
    let map = DofMap::build(ndof, &prescribed, &mpcs)?;
    let nfree = map.n_ind;
    if nfree == 0 {
        return err("RIKS: keine freien DOF.");
    }
    let (_, f_red) = map.reduce_inc(&[], &f_ext);
    let hist: Vec<Vec<plastic::GpHist>> = model
        .elements
        .iter()
        .map(|el| vec![plastic::GpHist::default(); plastic::n_gauss(el.kind).max(1)])
        .collect();

    let mut u_full = map.u0.clone();
    let mut lam = 0.0;
    let mut u_conv = u_full.clone();
    let mut lam_conv = 0.0;
    let mut du_prev = vec![0.0; nfree];
    let mut dl = 0.0;
    let mut total_iters = 0usize;
    let mut ninc = 0usize;
    let mut residual = 0.0;
    let mut solver_name = String::new();
    let mut last = assemble_nl(&model, &u_full, &hist)?;
    let mut retries = 0usize;
    let mut inc = 0usize;

    while inc < ctrl.max_inc.max(1) {
        last = assemble_nl(&model, &u_full, &hist)?;
        let (ff, _) = map.reduce_inc(&last.trips, &last.f_int);
        let vsol = match solve_kff(nfree, ff, &f_red) {
            Ok(s) => s,
            Err(_) => {
                retries += 1;
                if retries > 8 {
                    return err(format!(
                        "Riks: Tangente singulär bei λ={lam:.4} (Inkrement {inc})."
                    ));
                }
                dl *= 0.5;
                u_full.clone_from(&u_conv);
                lam = lam_conv;
                continue;
            }
        };
        solver_name = vsol.name.clone();
        let v = vsol.x;
        let vnorm_v = vnorm(&v).max(1e-30);
        let s = if inc == 0 || vdot(&v, &du_prev) >= 0.0 {
            1.0
        } else {
            -1.0
        };

        let mut dlam;
        let du_pred: Vec<f64>;
        if inc == 0 {
            dlam = ctrl.dlam.abs().clamp(ctrl.dlam_min, ctrl.dlam_max.max(ctrl.dlam_min));
            if lam + dlam > ctrl.period {
                dlam = (ctrl.period - lam).max(ctrl.dlam_min);
            }
            du_pred = v.iter().map(|x| dlam * x).collect();
            dl = vnorm(&du_pred).max(1e-18);
        } else {
            dlam = s * dl / vnorm_v;
            if dlam.abs() > ctrl.dlam_max {
                dlam = s * ctrl.dlam_max;
                dl = dlam.abs() * vnorm_v;
            }
            if dlam.abs() < ctrl.dlam_min {
                dlam = s * ctrl.dlam_min;
                dl = dlam.abs() * vnorm_v;
            }
            if dlam > 0.0 && lam + dlam > ctrl.period {
                dlam = ctrl.period - lam;
                du_pred = v.iter().map(|x| dlam * x).collect();
                dl = vnorm(&du_pred).max(1e-18);
            } else {
                du_pred = v.iter().map(|x| dlam * x).collect();
            }
        }

        let du_full = expand_red(&map, &du_pred, ndof);
        for i in 0..ndof {
            u_full[i] += du_full[i];
        }
        lam += dlam;
        let mut du_acc = du_pred;

        let mut conv = false;
        let mut it_count = 0usize;
        for it in 0..model.max_newton.max(1) {
            it_count = it + 1;
            total_iters += 1;
            last = assemble_nl(&model, &u_full, &hist)?;
            let mut r_full = vec![0.0; ndof];
            for i in 0..ndof {
                r_full[i] = lam * f_ext[i] - last.f_int[i];
            }
            let (ff, r_red) = map.reduce_inc(&last.trips, &r_full);
            residual = vnorm(&r_red);
            let fscale = (lam.abs() * vnorm(&f_red)).max(1.0);
            if residual < model.newton_tol * fscale {
                conv = true;
                break;
            }
            let du_i = match solve_kff(nfree, ff.clone(), &f_red) {
                Ok(s) => s,
                Err(_) => break,
            };
            let du_ii = match solve_kff(nfree, ff, &r_red) {
                Ok(s) => s,
                Err(_) => break,
            };
            solver_name = du_ii.name.clone();
            let dlam_c = crisfield_dlam(&du_i.x, &du_ii.x, &du_acc, dl);
            let mut du_c = vec![0.0; nfree];
            for i in 0..nfree {
                du_c[i] = du_ii.x[i] + dlam_c * du_i.x[i];
            }
            let du_c_full = expand_red(&map, &du_c, ndof);
            for i in 0..ndof {
                u_full[i] += du_c_full[i];
            }
            lam += dlam_c;
            for i in 0..nfree {
                du_acc[i] += du_c[i];
            }
            if vnorm(&du_c) < 1e-14 {
                conv = true;
                break;
            }
        }
        if !conv {
            retries += 1;
            if retries > 8 {
                return err(format!(
                    "Riks Inkrement {} konvergierte nicht (λ={lam:.4}, r={residual:.3e}).",
                    inc + 1
                ));
            }
            dl *= 0.5;
            u_full.clone_from(&u_conv);
            lam = lam_conv;
            continue;
        }
        retries = 0;
        inc += 1;
        ninc = inc;
        du_prev = du_acc;
        u_conv.clone_from(&u_full);
        lam_conv = lam;
        let n_des = 4.0;
        dl *= (n_des / (it_count as f64).max(1.0)).sqrt();
        if lam >= ctrl.period - 1e-8 {
            break;
        }
    }
    if ninc == 0 {
        return err("Riks: kein Inkrement konvergiert.");
    }

    last = assemble_nl(&model, &u_full, &hist)?;
    let mut r_full = vec![0.0; ndof];
    for i in 0..ndof {
        r_full[i] = lam * f_ext[i] - last.f_int[i];
    }
    let (_, r_red) = map.reduce_inc(&[], &r_full);
    residual = vnorm(&r_red);

    let solver = format!(
        "Riks (λ={lam:.4}, {ninc} inc, {solver_name}, {total_iters} iters)"
    );

    let mut u = vec![[0.0; 3]; nnode];
    let ur = vec![[0.0; 3]; nnode];
    let mut rf = vec![[0.0; 3]; nnode];
    let rm = vec![[0.0; 3]; nnode];
    for ni in 0..nnode {
        for d in 0..3 {
            u[ni][d] = u_full[dof_of(ndn, ni, d)];
            rf[ni][d] = last.f_int[dof_of(ndn, ni, d)] - lam * f_ext[dof_of(ndn, ni, d)];
        }
    }
    let mut accs = vec![[0.0; 6]; nnode];
    let mut cnt = vec![0.0; nnode];
    let mut gacc = vec![[0.0; 6]; nnode];
    let mut pacc = vec![0.0; nnode];
    let mut stress_gp = Vec::new();
    for (ei, el) in model.elements.iter().enumerate() {
        let s = last.cauchy[ei];
        let g = last.gl[ei];
        let p = last.peeq_e[ei];
        for &id in &el.nodes {
            let ni = model.node_index(id)?;
            for c in 0..6 {
                accs[ni][c] += s[c];
                gacc[ni][c] += g[c];
            }
            pacc[ni] += p;
            cnt[ni] += 1.0;
        }
        stress_gp.push((el.id, 1usize, s));
    }
    let mut stress = vec![[0.0; 6]; nnode];
    let mut strain = vec![[0.0; 6]; nnode];
    let mut vm = vec![0.0; nnode];
    let mut peeq = vec![0.0; nnode];
    for i in 0..nnode {
        if cnt[i] > 0.0 {
            for c in 0..6 {
                stress[i][c] = accs[i][c] / cnt[i];
                strain[i][c] = gacc[i][c] / cnt[i];
            }
            peeq[i] = pacc[i] / cnt[i];
        }
        vm[i] = von_mises(&stress[i]);
    }
    let frd_s = frd::write_frd(&model, &u, &stress, &rf, &strain, &peeq);
    let dat_s = dat::write_dat(&model, &u, &stress_gp, &rf);
    let procedure = model.procedure.name().to_string();
    let nsteps = model.steps.len().max(1);
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
        iters: total_iters,
        residual,
        time_ms: now_ms() - t0,
        procedure,
        frequencies: vec![],
        buckles: vec![],
        nsteps,
        lambda: lam,
        ninc,
        peeq,
    })
}
