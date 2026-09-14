//! Multi-point constraints: *EQUATION, *TIE, *RIGID BODY, *COUPLING.

use std::collections::{HashMap, HashSet};

use crate::error::{err, Result};
use crate::model::Model;

#[derive(Clone, Debug)]
pub struct Mpc {
    /// First term is the slave. Remaining terms: u_slave = sum (-c_i/c_s) u_i  (homogeneous).
    pub slave: usize,
    pub masters: Vec<(usize, f64)>, // (global dof, alpha)  u_s = sum alpha * u_m + u0_extra
    pub u0: f64,
}

#[derive(Clone, Debug)]
pub struct DofMap {
    pub n_full: usize,
    pub n_ind: usize,
    pub u0: Vec<f64>,
    pub t_row: Vec<Vec<(usize, f64)>>,
    pub ind_of: Vec<isize>,
}

impl DofMap {
    pub fn build(ndof: usize, prescribed: &HashMap<usize, f64>, mpcs: &[Mpc]) -> Result<Self> {
        let mut u0 = vec![0.0; ndof];
        for (&d, &v) in prescribed {
            if d < ndof {
                u0[d] = v;
            }
        }
        let mut slave_of: HashMap<usize, usize> = HashMap::new();
        for (i, m) in mpcs.iter().enumerate() {
            if m.slave >= ndof {
                continue;
            }
            if prescribed.contains_key(&m.slave) {
                continue;
            }
            if slave_of.insert(m.slave, i).is_some() {
                continue;
            }
        }

        // Expand nested slaves into independent / prescribed DOFs.
        let mut alphas: HashMap<usize, Vec<(usize, f64)>> = HashMap::new();
        let mut slave_u0: HashMap<usize, f64> = HashMap::new();
        for (si, mpc) in mpcs.iter().enumerate() {
            if !slave_of.get(&mpc.slave).map(|&k| k == si).unwrap_or(false) {
                continue;
            }
            let mut terms = mpc.masters.clone();
            let mut particular = mpc.u0;
            for _pass in 0..32 {
                let mut next = Vec::new();
                let mut changed = false;
                for (d, a) in terms {
                    if let Some(&k) = slave_of.get(&d) {
                        changed = true;
                        let nested = &mpcs[k];
                        particular += a * nested.u0;
                        for &(md, ma) in &nested.masters {
                            next.push((md, a * ma));
                        }
                    } else {
                        next.push((d, a));
                    }
                }
                terms = coalesce(next);
                if !changed {
                    break;
                }
            }
            // Split prescribed masters into particular solution.
            let mut ind_terms = Vec::new();
            for (d, a) in terms {
                if let Some(&v) = prescribed.get(&d) {
                    particular += a * v;
                } else if slave_of.contains_key(&d) {
                    return err("MPC: zyklische Mehrpunktbedingung.");
                } else {
                    ind_terms.push((d, a));
                }
            }
            alphas.insert(mpc.slave, ind_terms);
            slave_u0.insert(mpc.slave, particular);
            u0[mpc.slave] = particular;
        }

        let mut is_ind = vec![true; ndof];
        for &d in prescribed.keys() {
            if d < ndof {
                is_ind[d] = false;
            }
        }
        for &s in slave_of.keys() {
            if s < ndof {
                is_ind[s] = false;
            }
        }

        let mut ind_of = vec![-1isize; ndof];
        let mut n_ind = 0usize;
        for d in 0..ndof {
            if is_ind[d] {
                ind_of[d] = n_ind as isize;
                n_ind += 1;
            }
        }
        if n_ind == 0 && slave_of.is_empty() && prescribed.is_empty() {
            return err("Modell ohne freie Freiheitsgrade.");
        }

        let mut t_row = vec![Vec::new(); ndof];
        for d in 0..ndof {
            if is_ind[d] {
                t_row[d].push((ind_of[d] as usize, 1.0));
            } else if let Some(terms) = alphas.get(&d) {
                for &(md, a) in terms {
                    let col = ind_of[md];
                    if col >= 0 && a.abs() > 0.0 {
                        t_row[d].push((col as usize, a));
                    }
                }
            }
        }
        let _ = slave_u0;
        Ok(Self {
            n_full: ndof,
            n_ind,
            u0,
            t_row,
            ind_of,
        })
    }

    pub fn reduce(
        &self,
        trips: &[(usize, usize, f64)],
        f_full: &[f64],
    ) -> (Vec<(usize, usize, f64)>, Vec<f64>) {
        let mut rhs = vec![0.0; self.n_ind];
        for i in 0..self.n_full {
            if i >= f_full.len() {
                break;
            }
            let fi = f_full[i];
            if fi.abs() == 0.0 {
                continue;
            }
            for &(a, ta) in &self.t_row[i] {
                rhs[a] += ta * fi;
            }
        }
        let mut ff = Vec::new();
        for &(i, j, v) in trips {
            if i >= self.n_full || j >= self.n_full || v.abs() == 0.0 {
                continue;
            }
            let u0j = self.u0[j];
            for &(a, ta) in &self.t_row[i] {
                if u0j.abs() > 0.0 {
                    rhs[a] -= v * ta * u0j;
                }
                for &(b, tb) in &self.t_row[j] {
                    ff.push((a, b, v * ta * tb));
                }
            }
        }
        (ff, rhs)
    }

    /// Like [`reduce`], but for a Newton increment: `du = T dx` (no `u0` particular solution).
    pub fn reduce_inc(
        &self,
        trips: &[(usize, usize, f64)],
        r_full: &[f64],
    ) -> (Vec<(usize, usize, f64)>, Vec<f64>) {
        let mut rhs = vec![0.0; self.n_ind];
        for i in 0..self.n_full {
            if i >= r_full.len() {
                break;
            }
            let fi = r_full[i];
            if fi.abs() == 0.0 {
                continue;
            }
            for &(a, ta) in &self.t_row[i] {
                rhs[a] += ta * fi;
            }
        }
        let mut ff = Vec::new();
        for &(i, j, v) in trips {
            if i >= self.n_full || j >= self.n_full || v.abs() == 0.0 {
                continue;
            }
            for &(a, ta) in &self.t_row[i] {
                for &(b, tb) in &self.t_row[j] {
                    ff.push((a, b, v * ta * tb));
                }
            }
        }
        (ff, rhs)
    }

    pub fn reconstruct(&self, x: &[f64]) -> Vec<f64> {
        let mut u = self.u0.clone();
        for i in 0..self.n_full {
            for &(j, t) in &self.t_row[i] {
                if j < x.len() {
                    u[i] += t * x[j];
                }
            }
        }
        u
    }
}

fn coalesce(mut terms: Vec<(usize, f64)>) -> Vec<(usize, f64)> {
    terms.sort_by_key(|t| t.0);
    let mut out: Vec<(usize, f64)> = Vec::new();
    for (d, a) in terms {
        if let Some(last) = out.last_mut() {
            if last.0 == d {
                last.1 += a;
                continue;
            }
        }
        if a.abs() > 0.0 {
            out.push((d, a));
        }
    }
    out.retain(|t| t.1.abs() > 1e-16);
    out
}

pub fn equations_to_mpcs(model: &Model, ndn: usize) -> Result<Vec<Mpc>> {
    let mut mpcs = Vec::new();
    for eq in &model.equations {
        if eq.terms.len() < 2 {
            model_warn_skip();
            continue;
        }
        let mut terms: Vec<(usize, f64)> = Vec::new();
        for &(node, dof, coef) in &eq.terms {
            if dof >= ndn {
                continue;
            }
            let ni = model.node_index(node)?;
            terms.push((ndn * ni + dof, coef));
        }
        terms = coalesce(terms);
        if terms.len() < 2 {
            continue;
        }
        // Slave = largest |coef|
        let mut sidx = 0usize;
        for i in 1..terms.len() {
            if terms[i].1.abs() > terms[sidx].1.abs() {
                sidx = i;
            }
        }
        let (slave, cs) = terms[sidx];
        if cs.abs() < 1e-18 {
            return err("*EQUATION: Koeffizient des Slave-DOF ist null.");
        }
        let mut masters = Vec::new();
        for (i, (d, c)) in terms.iter().enumerate() {
            if i == sidx {
                continue;
            }
            masters.push((*d, -c / cs));
        }
        mpcs.push(Mpc {
            slave,
            masters,
            u0: eq.rhs / cs,
        });
    }
    Ok(mpcs)
}

fn model_warn_skip() {}

pub fn ties_to_mpcs(model: &Model, ndn: usize) -> Result<Vec<Mpc>> {
    let mut mpcs = Vec::new();
    let diag = bbox_diag(model);
    for tie in &model.ties {
        let slave_nodes = surface_nodes(model, &tie.slave)?;
        let master_nodes = surface_nodes(model, &tie.master)?;
        if slave_nodes.is_empty() || master_nodes.is_empty() {
            continue;
        }
        let tol = if tie.position_tol > 0.0 {
            tie.position_tol
        } else {
            1e-4 * diag.max(1.0)
        };
        let mut used_master = HashSet::new();
        for &s in &slave_nodes {
            let si = model.node_index(s)?;
            let xs = model.coords[si];
            let mut best = None;
            let mut best_d = f64::MAX;
            for &m in &master_nodes {
                if s == m {
                    continue;
                }
                let mi = model.node_index(m)?;
                let xm = model.coords[mi];
                let d = dist(xs, xm);
                if d < best_d {
                    best_d = d;
                    best = Some(m);
                }
            }
            let m = match best {
                Some(m) if best_d <= tol.max(1e-12) => m,
                Some(m) if !used_master.contains(&m) => m, // nearest even if far
                _ => continue,
            };
            used_master.insert(m);
            let mi = model.node_index(m)?;
            for dof in 0..3.min(ndn) {
                let slave = ndn * si + dof;
                let master = ndn * mi + dof;
                mpcs.push(Mpc {
                    slave,
                    masters: vec![(master, 1.0)],
                    u0: 0.0,
                });
            }
        }
    }
    Ok(mpcs)
}

pub fn rigid_to_mpcs(model: &Model, ndn: usize) -> Result<Vec<Mpc>> {
    let mut mpcs = Vec::new();
    if ndn < 6 && !model.rigid_bodies.is_empty() {
        return err("*RIGID BODY benötigt 6 DOF (Rotationen am Referenzknoten).");
    }
    for rb in &model.rigid_bodies {
        let refn = rb.ref_node;
        let ri = model.node_index(refn)?;
        let xr = model.coords[ri];
        let slaves = model.expand_nset(&rb.nset)?;
        for s in slaves {
            if s == refn {
                continue;
            }
            let si = model.node_index(s)?;
            let xs = model.coords[si];
            let rx = xs[0] - xr[0];
            let ry = xs[1] - xr[1];
            let rz = xs[2] - xr[2];
            // u_s = u_r + θ × r
            // ux: urx + θy*rz - θz*ry
            let urx = ndn * ri;
            let ury = ndn * ri + 1;
            let urz = ndn * ri + 2;
            let thx = ndn * ri + 3;
            let thy = ndn * ri + 4;
            let thz = ndn * ri + 5;
            mpcs.push(Mpc {
                slave: ndn * si,
                masters: vec![(urx, 1.0), (thy, rz), (thz, -ry)],
                u0: 0.0,
            });
            mpcs.push(Mpc {
                slave: ndn * si + 1,
                masters: vec![(ury, 1.0), (thz, rx), (thx, -rz)],
                u0: 0.0,
            });
            mpcs.push(Mpc {
                slave: ndn * si + 2,
                masters: vec![(urz, 1.0), (thx, ry), (thy, -rx)],
                u0: 0.0,
            });
            if ndn >= 6 {
                for r in 0..3 {
                    mpcs.push(Mpc {
                        slave: ndn * si + 3 + r,
                        masters: vec![(ndn * ri + 3 + r, 1.0)],
                        u0: 0.0,
                    });
                }
            }
        }
    }
    Ok(mpcs)
}

pub fn coupling_to_mpcs(model: &Model, ndn: usize) -> Result<Vec<Mpc>> {
    let mut mpcs = Vec::new();
    for c in &model.couplings {
        let nodes = surface_nodes(model, &c.surface)?;
        if nodes.is_empty() {
            continue;
        }
        let ri = model.node_index(c.ref_node)?;
        if c.kinematic {
            let xr = model.coords[ri];
            if ndn < 6 {
                return err("*COUPLING, KINEMATIC benötigt 6 DOF.");
            }
            for s in nodes {
                if s == c.ref_node {
                    continue;
                }
                let si = model.node_index(s)?;
                let xs = model.coords[si];
                let rx = xs[0] - xr[0];
                let ry = xs[1] - xr[1];
                let rz = xs[2] - xr[2];
                let urx = ndn * ri;
                let ury = ndn * ri + 1;
                let urz = ndn * ri + 2;
                let thx = ndn * ri + 3;
                let thy = ndn * ri + 4;
                let thz = ndn * ri + 5;
                if c.dofs.iter().any(|&d| d == 0) {
                    mpcs.push(Mpc {
                        slave: ndn * si,
                        masters: vec![(urx, 1.0), (thy, rz), (thz, -ry)],
                        u0: 0.0,
                    });
                }
                if c.dofs.iter().any(|&d| d == 1) {
                    mpcs.push(Mpc {
                        slave: ndn * si + 1,
                        masters: vec![(ury, 1.0), (thz, rx), (thx, -rz)],
                        u0: 0.0,
                    });
                }
                if c.dofs.iter().any(|&d| d == 2) {
                    mpcs.push(Mpc {
                        slave: ndn * si + 2,
                        masters: vec![(urz, 1.0), (thx, ry), (thy, -rx)],
                        u0: 0.0,
                    });
                }
            }
        } else {
            // Distributing: u_ref = (1/n) sum u_i  →  u_ref - (1/n) sum u_i = 0
            // Slave = ref node (so the independent DOFs stay on the surface).
            let n = nodes.len() as f64;
            if n < 1.0 {
                continue;
            }
            for &dof in &c.dofs {
                if dof >= ndn.min(3) {
                    continue;
                }
                let mut masters = Vec::new();
                for &s in &nodes {
                    if s == c.ref_node {
                        continue;
                    }
                    let si = model.node_index(s)?;
                    masters.push((ndn * si + dof, 1.0 / n));
                }
                mpcs.push(Mpc {
                    slave: ndn * ri + dof,
                    masters,
                    u0: 0.0,
                });
            }
        }
    }
    Ok(mpcs)
}

pub fn build_all_mpcs(model: &Model, ndn: usize) -> Result<Vec<Mpc>> {
    let mut v = equations_to_mpcs(model, ndn)?;
    v.extend(ties_to_mpcs(model, ndn)?);
    v.extend(rigid_to_mpcs(model, ndn)?);
    v.extend(coupling_to_mpcs(model, ndn)?);
    Ok(v)
}

fn surface_nodes(model: &Model, name: &str) -> Result<Vec<i32>> {
    let key = name.to_ascii_uppercase();
    if let Some(s) = model.surfaces.get(&key) {
        if !s.nodes.is_empty() {
            return Ok(s.nodes.clone());
        }
        let mut nodes = Vec::new();
        for &(eid, face) in &s.faces {
            if let Some(el) = model.elements.iter().find(|e| e.id == eid) {
                nodes.extend(face_nodes(el, face));
            }
        }
        nodes.sort();
        nodes.dedup();
        if !nodes.is_empty() {
            return Ok(nodes);
        }
    }
    model.expand_nset(&key)
}

pub(crate) fn face_nodes(el: &crate::model::Element, face: i32) -> Vec<i32> {
    use crate::model::ElemKind::*;
    let n = &el.nodes;
    let pick = |idx: &[usize]| -> Vec<i32> {
        idx.iter().filter_map(|&i| n.get(i).copied()).collect()
    };
    match (el.kind, face) {
        (Hex8 | Hex8I | Hex8R, 1) => pick(&[0, 1, 2, 3]),
        (Hex8 | Hex8I | Hex8R, 2) => pick(&[4, 7, 6, 5]),
        (Hex8 | Hex8I | Hex8R, 3) => pick(&[0, 4, 5, 1]),
        (Hex8 | Hex8I | Hex8R, 4) => pick(&[1, 5, 6, 2]),
        (Hex8 | Hex8I | Hex8R, 5) => pick(&[2, 6, 7, 3]),
        (Hex8 | Hex8I | Hex8R, 6) => pick(&[3, 7, 4, 0]),
        (Hex20 | Hex20R, f) if (1..=6).contains(&f) => {
            pick(&crate::quadratic::HEX20_FACE[(f - 1) as usize])
        }
        (Tet4, 1) => pick(&[0, 1, 2]),
        (Tet4, 2) => pick(&[0, 3, 1]),
        (Tet4, 3) => pick(&[1, 3, 2]),
        (Tet4, 4) => pick(&[2, 3, 0]),
        // C3D10 / C3D10T: 3 corners + 3 edge midsides (ccx faces).
        (Tet10 | Tet10T, 1) => pick(&[0, 1, 2, 4, 5, 6]),
        (Tet10 | Tet10T, 2) => pick(&[0, 3, 1, 7, 8, 4]),
        (Tet10 | Tet10T, 3) => pick(&[1, 3, 2, 8, 9, 5]),
        (Tet10 | Tet10T, 4) => pick(&[2, 3, 0, 9, 7, 6]),
        (Wedge6, 1) => pick(&[0, 1, 2]),
        (Wedge6, 2) => pick(&[3, 5, 4]),
        (Wedge6, 3) => pick(&[0, 1, 4, 3]),
        (Wedge6, 4) => pick(&[1, 2, 5, 4]),
        (Wedge6, 5) => pick(&[2, 0, 3, 5]),
        // C3D15: same connectivity as extra::wedge15_face_pressure.
        (Wedge15, 1) => pick(&[0, 1, 2, 6, 7, 8]),
        (Wedge15, 2) => pick(&[3, 5, 4, 11, 10, 9]),
        (Wedge15, 3) => pick(&[0, 1, 4, 3, 6, 13, 9, 12]),
        (Wedge15, 4) => pick(&[1, 2, 5, 4, 7, 14, 10, 13]),
        (Wedge15, 5) => pick(&[2, 0, 3, 5, 8, 12, 11, 14]),
        // Quadratic 2-D / shell edges S3–S6 (S1/S2 = the face = all nodes).
        (
            Quad8Ps | Quad8Pe | Quad8RPs | Quad8RPe | Cax8 | Cax8R | Shell8 | Shell8R | Mem8,
            3,
        ) => pick(&[0, 1, 4]),
        (
            Quad8Ps | Quad8Pe | Quad8RPs | Quad8RPe | Cax8 | Cax8R | Shell8 | Shell8R | Mem8,
            4,
        ) => pick(&[1, 2, 5]),
        (
            Quad8Ps | Quad8Pe | Quad8RPs | Quad8RPe | Cax8 | Cax8R | Shell8 | Shell8R | Mem8,
            5,
        ) => pick(&[2, 3, 6]),
        (
            Quad8Ps | Quad8Pe | Quad8RPs | Quad8RPe | Cax8 | Cax8R | Shell8 | Shell8R | Mem8,
            6,
        ) => pick(&[3, 0, 7]),
        (Tri6Ps | Tri6Pe | Cax6, 1) => pick(&[0, 1, 3]),
        (Tri6Ps | Tri6Pe | Cax6, 2) => pick(&[1, 2, 4]),
        (Tri6Ps | Tri6Pe | Cax6, 3) => pick(&[2, 0, 5]),
        _ => n.clone(),
    }
}

fn bbox_diag(model: &Model) -> f64 {
    if model.coords.is_empty() {
        return 1.0;
    }
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
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Rotate element stiffness into local nodal axes: ke_local = T^T ke_global T.
pub fn transform_ke(
    ke: &mut [f64],
    nn: usize,
    local_dim: usize,
    nodes: &[i32],
    transforms: &HashMap<i32, [[f64; 3]; 3]>,
) {
    if transforms.is_empty() {
        return;
    }
    let nd = nn * local_dim;
    if ke.len() < nd * nd {
        return;
    }
    // T is block-diagonal of 3x3 (and another 3x3 for rotations if local_dim>=6)
    let mut t = vec![0.0; nd * nd];
    for a in 0..nn {
        let r = transforms.get(&nodes[a]).copied().unwrap_or([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        // columns of R are local axes in global: u_g = R u_l
        for p in 0..3.min(local_dim) {
            for q in 0..3.min(local_dim) {
                t[(a * local_dim + p) * nd + (a * local_dim + q)] = r[p][q];
            }
        }
        if local_dim >= 6 {
            for p in 0..3 {
                for q in 0..3 {
                    t[(a * local_dim + 3 + p) * nd + (a * local_dim + 3 + q)] = r[p][q];
                }
            }
        }
    }
    // ke := T^T ke T
    let mut tmp = vec![0.0; nd * nd];
    for i in 0..nd {
        for j in 0..nd {
            let mut s = 0.0;
            for k in 0..nd {
                s += ke[i * nd + k] * t[k * nd + j];
            }
            tmp[i * nd + j] = s;
        }
    }
    for i in 0..nd {
        for j in 0..nd {
            let mut s = 0.0;
            for k in 0..nd {
                s += t[k * nd + i] * tmp[k * nd + j];
            }
            ke[i * nd + j] = s;
        }
    }
}

pub fn dofs_to_global(
    u: &mut [f64],
    ndn: usize,
    node_ids: &[i32],
    transforms: &HashMap<i32, [[f64; 3]; 3]>,
) {
    if transforms.is_empty() {
        return;
    }
    for (ni, id) in node_ids.iter().enumerate() {
        let Some(r) = transforms.get(id) else {
            continue;
        };
        let base = ndn * ni;
        if base + 2 >= u.len() {
            continue;
        }
        let ul = [u[base], u[base + 1], u[base + 2]];
        for p in 0..3 {
            u[base + p] = r[p][0] * ul[0] + r[p][1] * ul[1] + r[p][2] * ul[2];
        }
        if ndn >= 6 && base + 5 < u.len() {
            let rl = [u[base + 3], u[base + 4], u[base + 5]];
            for p in 0..3 {
                u[base + 3 + p] = r[p][0] * rl[0] + r[p][1] * rl[1] + r[p][2] * rl[2];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::face_nodes;
    use crate::model::{ElemKind, Element};

    fn el(kind: ElemKind, n: usize) -> Element {
        Element {
            id: 1,
            kind,
            nodes: (1..=n as i32).collect(),
            elset: String::new(),
        }
    }

    #[test]
    fn c3d20_faces_include_midsides() {
        let e = el(ElemKind::Hex20, 20);
        assert_eq!(face_nodes(&e, 1), vec![1, 2, 3, 4, 9, 10, 11, 12]);
        assert_eq!(face_nodes(&e, 2), vec![5, 8, 7, 6, 16, 15, 14, 13]);
        assert_eq!(face_nodes(&e, 3).len(), 8);
        assert_eq!(face_nodes(&el(ElemKind::Hex20R, 20), 1), face_nodes(&e, 1));
    }

    #[test]
    fn c3d10_faces_include_midsides() {
        let e = el(ElemKind::Tet10, 10);
        assert_eq!(face_nodes(&e, 1), vec![1, 2, 3, 5, 6, 7]);
        assert_eq!(face_nodes(&e, 2), vec![1, 4, 2, 8, 9, 5]);
        assert_eq!(face_nodes(&el(ElemKind::Tet10T, 10), 1), face_nodes(&e, 1));
        // linear tet stays 3-node
        assert_eq!(face_nodes(&el(ElemKind::Tet4, 4), 1), vec![1, 2, 3]);
    }

    #[test]
    fn c3d15_faces_include_midsides() {
        let e = el(ElemKind::Wedge15, 15);
        assert_eq!(face_nodes(&e, 1), vec![1, 2, 3, 7, 8, 9]);
        assert_eq!(face_nodes(&e, 3).len(), 8);
        assert_eq!(face_nodes(&el(ElemKind::Wedge6, 6), 1), vec![1, 2, 3]);
        assert_eq!(face_nodes(&el(ElemKind::Wedge6, 6), 3), vec![1, 2, 5, 4]);
    }

    #[test]
    fn shell8_spos_keeps_all_eight() {
        let e = el(ElemKind::Shell8, 8);
        assert_eq!(face_nodes(&e, 1).len(), 8, "SPOS/S1 is the 8-node face");
        assert_eq!(face_nodes(&e, 3), vec![1, 2, 5]);
    }
}
