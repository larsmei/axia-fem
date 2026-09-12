//! Penalty node-to-surface contact (`*CONTACT PAIR`, frictionless).
//! Gap g = n·(x_s − x_c); active if g < 0. Force on slave F_s = −k n g.

use crate::constraint;
use crate::error::{err, Result};
use crate::model::{ContactPair, Model};

pub struct ContactForce {
    pub trips: Vec<(usize, usize, f64)>,
    pub f: Vec<f64>,
    pub n_active: usize,
}

fn vsub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn vadd(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn vscale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn vdot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn vcross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn vnorm(a: [f64; 3]) -> f64 {
    vdot(a, a).sqrt()
}

fn pos(model: &Model, u: &[f64], ndn: usize, id: i32) -> Result<[f64; 3]> {
    let ni = model.node_index(id)?;
    let x = model.coords[ni];
    Ok([
        x[0] + u.get(ndn * ni).copied().unwrap_or(0.0),
        x[1] + u.get(ndn * ni + 1).copied().unwrap_or(0.0),
        x[2] + u.get(ndn * ni + 2).copied().unwrap_or(0.0),
    ])
}

fn dof_t(ndn: usize, ni: usize, d: usize) -> usize {
    ndn * ni + d
}

struct Hit {
    gap: f64,
    n: [f64; 3],
    master: Vec<i32>,
    shape: Vec<f64>,
}

fn project_tri(xs: [f64; 3], p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> Option<(f64, [f64; 3], [f64; 3])> {
    let e1 = vsub(p1, p0);
    let e2 = vsub(p2, p0);
    let mut n = vcross(e1, e2);
    let a2 = vnorm(n);
    if a2 < 1e-18 {
        return None;
    }
    n = vscale(n, 1.0 / a2);
    let gap = vdot(n, vsub(xs, p0));
    let q = vsub(xs, vscale(n, gap));
    let v0 = e1;
    let v1 = e2;
    let v2 = vsub(q, p0);
    let d00 = vdot(v0, v0);
    let d01 = vdot(v0, v1);
    let d11 = vdot(v1, v1);
    let d20 = vdot(v2, v0);
    let d21 = vdot(v2, v1);
    let den = d00 * d11 - d01 * d01;
    if den.abs() < 1e-30 {
        return None;
    }
    let v = (d11 * d20 - d01 * d21) / den;
    let w = (d00 * d21 - d01 * d20) / den;
    let u = 1.0 - v - w;
    const EPS: f64 = 1e-7;
    if u < -EPS || v < -EPS || w < -EPS {
        return None;
    }
    Some((gap, [u, v, w], n))
}

fn outward(n: [f64; 3], face_c: [f64; 3], elem_c: [f64; 3]) -> [f64; 3] {
    if vdot(n, vsub(face_c, elem_c)) < 0.0 {
        vscale(n, -1.0)
    } else {
        n
    }
}

fn elem_centroid(model: &Model, eid: i32) -> Option<[f64; 3]> {
    let el = model.elements.iter().find(|e| e.id == eid)?;
    let mut c = [0.0; 3];
    let mut n = 0.0;
    for &id in &el.nodes {
        if let Ok(i) = model.node_index(id) {
            let x = model.coords[i];
            c = vadd(c, x);
            n += 1.0;
        }
    }
    if n > 0.0 {
        Some(vscale(c, 1.0 / n))
    } else {
        None
    }
}

fn master_faces(model: &Model, name: &str) -> Result<Vec<(i32, Vec<i32>)>> {
    let key = name.to_ascii_uppercase();
    let s = model
        .surfaces
        .get(&key)
        .ok_or_else(|| crate::error::FemError(format!("*CONTACT: Oberfläche {name} fehlt.")))?;
    if s.faces.is_empty() {
        return err(format!(
            "*CONTACT: Master {name} braucht TYPE=ELEMENT (Flächen)."
        ));
    }
    let mut out = Vec::new();
    for &(eid, face) in &s.faces {
        let el = model
            .elements
            .iter()
            .find(|e| e.id == eid)
            .ok_or_else(|| crate::error::FemError(format!("*CONTACT: Element {eid} fehlt.")))?;
        let nodes = constraint::face_nodes(el, face);
        if nodes.len() >= 3 {
            out.push((eid, nodes));
        }
    }
    if out.is_empty() {
        return err(format!("*CONTACT: Master {name} hat keine gültigen Flächen."));
    }
    Ok(out)
}

fn slave_ids(model: &Model, name: &str) -> Result<Vec<i32>> {
    let key = name.to_ascii_uppercase();
    if let Some(s) = model.surfaces.get(&key) {
        let mut ids = s.nodes.clone();
        for &(eid, face) in &s.faces {
            if let Some(el) = model.elements.iter().find(|e| e.id == eid) {
                ids.extend(constraint::face_nodes(el, face));
            }
        }
        ids.sort();
        ids.dedup();
        if !ids.is_empty() {
            return Ok(ids);
        }
    }
    model.expand_nset(&key)
}

fn find_hit(model: &Model, u: &[f64], ndn: usize, xs: [f64; 3], faces: &[(i32, Vec<i32>)]) -> Result<Option<Hit>> {
    let mut best: Option<Hit> = None;
    for (eid, mnodes) in faces {
        let mut pts = Vec::with_capacity(mnodes.len());
        for &id in mnodes {
            pts.push(pos(model, u, ndn, id)?);
        }
        let npt = pts.len();
        let mut face_c = [0.0; 3];
        for p in &pts {
            face_c = vadd(face_c, *p);
        }
        face_c = vscale(face_c, 1.0 / npt as f64);
        let elem_c = elem_centroid(model, *eid).unwrap_or(face_c);

        let tris: Vec<[usize; 3]> = if npt == 3 {
            vec![[0, 1, 2]]
        } else if npt >= 4 {
            vec![[0, 1, 2], [0, 2, 3]]
        } else {
            continue;
        };
        for t in tris {
            if let Some((mut gap, bary, mut n)) = project_tri(xs, pts[t[0]], pts[t[1]], pts[t[2]]) {
                let n_out = outward(n, face_c, elem_c);
                if vdot(n_out, n) < 0.0 {
                    n = n_out;
                    gap = -gap;
                }
                let mut shape = vec![0.0; npt];
                shape[t[0]] = bary[0];
                shape[t[1]] = bary[1];
                shape[t[2]] = bary[2];
                let better = match &best {
                    None => true,
                    Some(h) => gap.abs() < h.gap.abs(),
                };
                if better {
                    best = Some(Hit {
                        gap,
                        n,
                        master: mnodes.clone(),
                        shape,
                    });
                }
            }
        }
    }
    Ok(best)
}

pub fn assemble(
    model: &Model,
    ndn: usize,
    ndof: usize,
    u: &[f64],
) -> Result<ContactForce> {
    let mut trips = Vec::new();
    let mut f = vec![0.0; ndof];
    let mut n_active = 0usize;
    for pair in &model.contact_pairs {
        add_pair(model, pair, ndn, u, &mut trips, &mut f, &mut n_active)?;
    }
    Ok(ContactForce {
        trips,
        f,
        n_active,
    })
}

fn add_pair(
    model: &Model,
    pair: &ContactPair,
    ndn: usize,
    u: &[f64],
    trips: &mut Vec<(usize, usize, f64)>,
    f: &mut [f64],
    n_active: &mut usize,
) -> Result<()> {
    let kn = pair.kn.max(0.0);
    if kn == 0.0 {
        return Ok(());
    }
    let faces = master_faces(model, &pair.master)?;
    let slaves = slave_ids(model, &pair.slave)?;
    let mut master_set = std::collections::HashSet::new();
    for (_, n) in &faces {
        for &id in n {
            master_set.insert(id);
        }
    }
    for sid in slaves {
        if master_set.contains(&sid) {
            continue;
        }
        let xs = pos(model, u, ndn, sid)?;
        let Some(hit) = find_hit(model, u, ndn, xs, &faces)? else {
            continue;
        };
        if hit.gap >= 0.0 {
            continue;
        }
        *n_active += 1;
        let sni = model.node_index(sid)?;
        let mut nodes = vec![sid];
        let mut w = vec![-1.0];
        for (a, &mid) in hit.master.iter().enumerate() {
            nodes.push(mid);
            w.push(hit.shape.get(a).copied().unwrap_or(0.0));
        }
        // F_s = −kn g n, F_m = kn g N n. Residual contribution is −F_int, but
        // here `f` is the internal contact force (same sign as F).
        let fnod = -kn * hit.gap;
        for (a, &id) in nodes.iter().enumerate() {
            let ni = model.node_index(id)?;
            for d in 0..3 {
                f[dof_t(ndn, ni, d)] += fnod * w[a] * hit.n[d];
            }
        }
        for (a, &ida) in nodes.iter().enumerate() {
            let ia = model.node_index(ida)?;
            for (b, &idb) in nodes.iter().enumerate() {
                let ib = model.node_index(idb)?;
                let kab = kn * w[a] * w[b];
                for i in 0..3 {
                    for j in 0..3 {
                        let v = kab * hit.n[i] * hit.n[j];
                        if v.abs() > 0.0 {
                            trips.push((dof_t(ndn, ia, i), dof_t(ndn, ib, j), v));
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn project_inside_unit_tri() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        let xs = [0.25, 0.25, 0.5];
        let (g, b, n) = project_tri(xs, p0, p1, p2).unwrap();
        assert!((g - 0.5).abs() < 1e-12, "gap={g}");
        assert!((n[2] - 1.0).abs() < 1e-12);
        assert!((b[0] - 0.5).abs() < 1e-12);
    }
}
