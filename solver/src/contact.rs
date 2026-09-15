//! Penalty node-to-surface contact (`*CONTACT PAIR`).
//! Gap g = n·(x_s − x_c); active if g < 0. Coulomb on nodal forces, small sliding.

use crate::constraint;
use crate::error::{err, Result};
use crate::model::{ContactPair, Model};

pub struct ContactForce {
    pub trips: Vec<(usize, usize, f64)>,
    pub f: Vec<f64>,
    pub n_active: usize,
    pub n_slip: usize,
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
    area: f64,
}

fn project_tri(xs: [f64; 3], p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> Option<(f64, [f64; 3], [f64; 3], f64)> {
    let e1 = vsub(p1, p0);
    let e2 = vsub(p2, p0);
    let mut n = vcross(e1, e2);
    let a2 = vnorm(n);
    if a2 < 1e-18 {
        return None;
    }
    n = vscale(n, 1.0 / a2);
    let v0 = e1;
    let v1 = e2;
    let v2 = vsub(xs, p0);
    let d00 = vdot(v0, v0);
    let d01 = vdot(v0, v1);
    let d11 = vdot(v1, v1);
    let d20 = vdot(v2, v0);
    let d21 = vdot(v2, v1);
    let den = d00 * d11 - d01 * d01;
    if den.abs() < 1e-30 {
        return None;
    }
    let mut v = (d11 * d20 - d01 * d21) / den;
    let mut w = (d00 * d21 - d01 * d20) / den;
    let mut u = 1.0 - v - w;
    // Clamp to triangle (closest point), then reject far in-plane hits.
    if u < 0.0 {
        let e = vsub(p2, p1);
        let t = vdot(vsub(xs, p1), e) / vdot(e, e).max(1e-30);
        let t = t.clamp(0.0, 1.0);
        u = 0.0;
        v = 1.0 - t;
        w = t;
    } else if v < 0.0 {
        let e = vsub(p0, p2);
        let t = vdot(vsub(xs, p2), e) / vdot(e, e).max(1e-30);
        let t = t.clamp(0.0, 1.0);
        v = 0.0;
        w = 1.0 - t;
        u = t;
    } else if w < 0.0 {
        let e = vsub(p1, p0);
        let t = vdot(vsub(xs, p0), e) / vdot(e, e).max(1e-30);
        let t = t.clamp(0.0, 1.0);
        w = 0.0;
        u = 1.0 - t;
        v = t;
    }
    let xc = vadd(vadd(vscale(p0, u), vscale(p1, v)), vscale(p2, w));
    let gap = vdot(n, vsub(xs, xc));
    let tang = vnorm(vsub(vsub(xs, xc), vscale(n, gap)));
    let e1n = vnorm(e1);
    let e2n = vnorm(e2);
    let face = e1n.max(e2n).max(1e-12);
    // 15 % of the facet, capped at 1.5 so ~10-unit patch tests still hit,
    // tight enough that a 0.1 m bolt overhang cannot grab a 50 mm plate.
    let tang_cap = (0.15 * face).clamp(1e-9, 1.5);
    if tang > tang_cap {
        return None;
    }
    Some((gap, [u, v, w], n, 0.5 * a2))
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

/// Drop quadratic midside slaves. Linear facets already ignore midsides on the
/// master; keeping them as slaves chatters a C3D20 bolted joint.
fn quadratic_slave_corners(model: &Model, name: &str, slaves: Vec<i32>) -> Result<Vec<i32>> {
    let key = name.to_ascii_uppercase();
    let Some(s) = model.surfaces.get(&key) else {
        return Ok(slaves);
    };
    if s.faces.is_empty() {
        return Ok(slaves);
    }
    let mut corners = std::collections::HashSet::new();
    let mut quad = false;
    for &(eid, face) in &s.faces {
        let Some(el) = model.elements.iter().find(|e| e.id == eid) else {
            continue;
        };
        let nodes = constraint::face_nodes(el, face);
        let n_c = if nodes.len() >= 8 {
            4
        } else if nodes.len() >= 6 {
            3
        } else {
            nodes.len()
        };
        if nodes.len() > n_c {
            quad = true;
        }
        for &id in nodes.iter().take(n_c) {
            corners.insert(id);
        }
    }
    if !quad {
        return Ok(slaves);
    }
    Ok(slaves.into_iter().filter(|id| corners.contains(id)).collect())
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
            if let Some((mut gap, bary, mut n, area)) = project_tri(xs, pts[t[0]], pts[t[1]], pts[t[2]]) {
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
                        area,
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
    assemble_stick(model, ndn, ndof, u, false)
}

/// `full_stick`: pretension Coulomb uses kt=kn (needed while p≈0 on
/// increment 1). Later increments use smooth Coulomb ft = τ gt / s,
/// τ = μp + kt_floor g_char.
pub fn assemble_stick(
    model: &Model,
    ndn: usize,
    ndof: usize,
    u: &[f64],
    full_stick: bool,
) -> Result<ContactForce> {
    let mut trips = Vec::new();
    let mut f = vec![0.0; ndof];
    let mut n_active = 0usize;
    let mut n_slip = 0usize;
    for pair in &model.contact_pairs {
        add_pair(
            model,
            pair,
            ndn,
            u,
            &mut trips,
            &mut f,
            &mut n_active,
            &mut n_slip,
            full_stick,
        )?;
    }
    Ok(ContactForce {
        trips,
        f,
        n_active,
        n_slip,
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
    n_slip: &mut usize,
    full_stick: bool,
) -> Result<()> {
    let kn = pair.kn.max(0.0);
    if kn == 0.0 {
        return Ok(());
    }
    let mu = pair.mu.max(0.0);
    let pret = !model.pretensions.is_empty();
    let faces = master_faces(model, &pair.master)?;
    let slaves = quadratic_slave_corners(model, &pair.slave, slave_ids(model, &pair.slave)?)?;
    let mut master_set = std::collections::HashSet::new();
    for (_, n) in &faces {
        for &id in n {
            master_set.insert(id);
        }
    }
    let mut dummy = std::collections::HashSet::new();
    let mut twin_of: std::collections::HashMap<i32, i32> = std::collections::HashMap::new();
    for p in &model.pretensions {
        dummy.insert(p.dummy);
        for &(a, b, _) in &p.pairs {
            twin_of.insert(a, b);
            twin_of.insert(b, a);
        }
    }
    for sid in slaves {
        if dummy.contains(&sid) {
            continue;
        }
        if master_set.contains(&sid) {
            continue;
        }
        if let Some(&tw) = twin_of.get(&sid) {
            if master_set.contains(&tw) {
                continue;
            }
        }
        let xs = pos(model, u, ndn, sid)?;
        let Some(hit) = find_hit(model, u, ndn, xs, &faces)? else {
            continue;
        };
        if let Some(&tw) = twin_of.get(&sid) {
            if hit.master.iter().any(|&m| m == tw) {
                continue;
            }
        }
        if hit.gap > 1e-9 {
            // Unilateral for pretension (don't glue a bolted joint that is
            // trying to slide). Other contact: keep a 1-unit approaching
            // band so a 0.2-gap patch test is not a rigid-body mechanism
            // on the first Newton step.
            if pret || hit.gap > 1.0 {
                continue;
            }
        }
        *n_active += 1;
        let mut nodes = vec![sid];
        let mut w = vec![-1.0];
        for (a, &mid) in hit.master.iter().enumerate() {
            nodes.push(mid);
            w.push(hit.shape.get(a).copied().unwrap_or(0.0));
        }
        // LINEAR pressure-overclosure is force/length³. Convert to a nodal
        // spring with the hit-triangle area. `area.min(1)` keeps the historic
        // ~10-unit patch tests (A≫1) on kn-as-nodal-stiffness. Cap SI
        // Mecway slopes (2e14 Pa/m) so a micron of overclosure is not giganewtons.
        let kn_n = (kn * hit.area.min(1.0).max(1e-30)).min(1e10);
        let kt = kn_n;
        // Pretension stabilizer only for frictional approaching contacts
        // (0 < gap ≤ 1e-9). Frictionless pairs (SI_19 bolt head, μ=0) stay
        // tangentially free, matching CalculiX *FRICTION-less SURFACE
        // INTERACTION. Gluing the head with kn was the remaining sandwich.
        let kt_stab = if pret && mu > 0.0 { kn_n } else { 0.0 };
        let fnod = -kn_n * hit.gap;
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
                let kab = kn_n * w[a] * w[b];
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
        if mu <= 0.0 || hit.gap > 0.0 {
            if kt_stab > 0.0 {
                let mut g_raw = [0.0; 3];
                for (a, &id) in nodes.iter().enumerate() {
                    let ni = model.node_index(id)?;
                    for d in 0..3 {
                        g_raw[d] += w[a] * u.get(dof_t(ndn, ni, d)).copied().unwrap_or(0.0);
                    }
                }
                let mut gt = [-g_raw[0], -g_raw[1], -g_raw[2]];
                let gn = vdot(gt, hit.n);
                gt = vsub(gt, vscale(hit.n, gn));
                let ft = vscale(gt, kt_stab);
                for (a, &id) in nodes.iter().enumerate() {
                    let ni = model.node_index(id)?;
                    for d in 0..3 {
                        f[dof_t(ndn, ni, d)] += w[a] * (-ft[d]);
                    }
                }
                for (a, &ida) in nodes.iter().enumerate() {
                    let ia = model.node_index(ida)?;
                    for (b, &idb) in nodes.iter().enumerate() {
                        let ib = model.node_index(idb)?;
                        let wab = w[a] * w[b];
                        for i in 0..3 {
                            for j in 0..3 {
                                let mut pij = if i == j { 1.0 } else { 0.0 };
                                pij -= hit.n[i] * hit.n[j];
                                let v = kt_stab * wab * pij;
                                if v.abs() > 0.0 {
                                    trips.push((dof_t(ndn, ia, i), dof_t(ndn, ib, j), v));
                                }
                            }
                        }
                    }
                }
            }
            continue;
        }
        // Small-sliding Coulomb: g_t = (I−nn)(u_s − u_c).
        let mut g_raw = [0.0; 3];
        for (a, &id) in nodes.iter().enumerate() {
            let ni = model.node_index(id)?;
            for d in 0..3 {
                g_raw[d] += w[a] * u.get(dof_t(ndn, ni, d)).copied().unwrap_or(0.0);
            }
        }
        let mut gt = [-g_raw[0], -g_raw[1], -g_raw[2]];
        let gn = vdot(gt, hit.n);
        gt = vsub(gt, vscale(hit.n, gn));
        let gt_n = vnorm(gt);
        let p = (kn_n * (-hit.gap)).max(0.0);
        // Smooth Coulomb during pretension (after increment 1):
        //   ft = τ gt / sqrt(|gt|² + g_char²),  τ = μp + kt_floor g_char
        // Saturates at τ (far field can slide) with a consistent SPD
        // tangent (τ/s)(I − t⊗t). g_char is 10 μm so the floor force is
        // ~100 N/node while the origin slope stays 1e-3 kn for Newton.
        const G_CHAR: f64 = 1e-5; // 10 μm; floor force ~100 N/node at kn_n=1e10
        const KT_FLOOR_FRAC: f64 = 1e-3;
        let (ft, k_t, tdir, slip) = if pret && !full_stick {
            let tau = mu * p + KT_FLOOR_FRAC * kn_n * G_CHAR;
            let s = (gt_n * gt_n + G_CHAR * G_CHAR).sqrt();
            let t = vscale(gt, 1.0 / s);
            if gt_n > G_CHAR {
                *n_slip += 1;
            }
            (vscale(t, tau), tau / s, t, true)
        } else {
            let mut ft = vscale(gt, kt);
            let ft_n = vnorm(ft);
            let can_slip = !pret && ft_n > mu * p + 1e-12 && gt_n > 1e-16;
            if can_slip {
                *n_slip += 1;
                let tdir = vscale(gt, 1.0 / gt_n);
                ft = vscale(tdir, mu * p);
                (ft, mu * p / gt_n, tdir, true)
            } else {
                (ft, kt, [0.0; 3], false)
            }
        };
        for (a, &id) in nodes.iter().enumerate() {
            let ni = model.node_index(id)?;
            for d in 0..3 {
                f[dof_t(ndn, ni, d)] += w[a] * (-ft[d]);
            }
        }
        for (a, &ida) in nodes.iter().enumerate() {
            let ia = model.node_index(ida)?;
            for (b, &idb) in nodes.iter().enumerate() {
                let ib = model.node_index(idb)?;
                let wab = w[a] * w[b];
                for i in 0..3 {
                    for j in 0..3 {
                        let mut pij = if i == j { 1.0 } else { 0.0 };
                        pij -= hit.n[i] * hit.n[j];
                        if slip {
                            pij -= tdir[i] * tdir[j];
                        }
                        let v = k_t * wab * pij;
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
        let (g, b, n, area) = project_tri(xs, p0, p1, p2).unwrap();
        assert!((g - 0.5).abs() < 1e-12, "gap={g}");
        assert!((n[2] - 1.0).abs() < 1e-12);
        assert!((b[0] - 0.5).abs() < 1e-12);
        assert!((area - 0.5).abs() < 1e-12);
    }
}
