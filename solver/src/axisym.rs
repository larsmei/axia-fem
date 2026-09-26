//! Axisymmetric solids CAX3/CAX4/CAX4R/CAX6/CAX8/CAX8R.
//! r = x, z = y. Strains [εr, εz, εθ, γrz], 2 DOF/node (ur, uz).
//! Integration weight 2π r det J (CalculiX convention).

use crate::elem::{gemm_bt_d_b, invert2, sigma_from_b, G2};
use crate::error::{err, Result};
use crate::quadratic::{quad8_dndx, QUAD8_XI, G3, W3};

const QUAD_XI: [[f64; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];

fn d_axisym(e: f64, nu: f64) -> Result<[f64; 16]> {
    if nu >= 0.5 || nu <= -1.0 {
        return err(format!("Ungültige Querkontraktion nu={nu}"));
    }
    if e <= 0.0 {
        return err(format!("Ungültiger E-Modul E={e}"));
    }
    let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let d11 = lam + 2.0 * mu;
    let mut d = [0.0; 16];
    for i in 0..3 {
        for j in 0..3 {
            d[i * 4 + j] = if i == j { d11 } else { lam };
        }
    }
    d[15] = mu;
    Ok(d)
}

fn fill_b_cax(b: &mut [f64], nnode: usize, dndx: &[[f64; 2]], nshp: &[f64], r: f64) {
    // B 4 × (2 nnode): εr, εz, εθ, γrz
    let n = 2 * nnode;
    let invr = if r.abs() > 1e-12 { 1.0 / r } else { 0.0 };
    for i in 0..nnode {
        let c = 2 * i;
        let dr = dndx[i][0];
        let dz = dndx[i][1];
        b[0 * n + c] = dr; // εr = dur/dr
        b[1 * n + c + 1] = dz; // εz = duz/dz
        b[2 * n + c] = if r.abs() > 1e-12 {
            nshp[i] * invr
        } else {
            dr
        }; // εθ = ur/r (L'Hôpital on axis)
        b[3 * n + c] = dz; // γrz
        b[3 * n + c + 1] = dr;
    }
}

fn quad4_shape(xi: f64, eta: f64) -> ([f64; 4], [[f64; 2]; 4]) {
    let mut n = [0.0; 4];
    let mut dn = [[0.0; 2]; 4];
    for i in 0..4 {
        let x = QUAD_XI[i][0];
        let e = QUAD_XI[i][1];
        n[i] = 0.25 * (1.0 + x * xi) * (1.0 + e * eta);
        dn[i][0] = 0.25 * x * (1.0 + e * eta);
        dn[i][1] = 0.25 * e * (1.0 + x * xi);
    }
    (n, dn)
}

fn quad4_dndx(xy: &[[f64; 2]], xi: f64, eta: f64) -> Result<([[f64; 2]; 4], f64, [f64; 4])> {
    let (n, dn) = quad4_shape(xi, eta);
    let mut j = [[0.0; 2]; 2];
    for a in 0..4 {
        j[0][0] += dn[a][0] * xy[a][0];
        j[0][1] += dn[a][1] * xy[a][0];
        j[1][0] += dn[a][0] * xy[a][1];
        j[1][1] += dn[a][1] * xy[a][1];
    }
    let (inv, det) = invert2(j)?;
    let mut dndx = [[0.0; 2]; 4];
    for a in 0..4 {
        dndx[a][0] = inv[0][0] * dn[a][0] + inv[1][0] * dn[a][1];
        dndx[a][1] = inv[0][1] * dn[a][0] + inv[1][1] * dn[a][1];
    }
    Ok((dndx, det, n))
}

fn radius(xy: &[[f64; 2]], nshp: &[f64], nn: usize) -> f64 {
    let mut r = 0.0;
    for a in 0..nn {
        r += nshp[a] * xy[a][0];
    }
    r.max(0.0)
}

fn two_pi() -> f64 {
    2.0 * std::f64::consts::PI
}

fn cax4_gauss(reduced: bool) -> Vec<(f64, f64, f64)> {
    if reduced {
        vec![(0.0, 0.0, 4.0)]
    } else {
        let mut o = Vec::new();
        for &xi in &[-G2, G2] {
            for &eta in &[-G2, G2] {
                o.push((xi, eta, 1.0));
            }
        }
        o
    }
}

fn cax8_gauss(reduced: bool) -> Vec<(f64, f64, f64)> {
    let mut o = Vec::new();
    if reduced {
        for &xi in &[-G2, G2] {
            for &eta in &[-G2, G2] {
                o.push((xi, eta, 1.0));
            }
        }
    } else {
        for i in 0..3 {
            for j in 0..3 {
                o.push((G3[i], G3[j], W3[i] * W3[j]));
            }
        }
    }
    o
}

pub fn cax4_stiffness(xyz: &[[f64; 3]], e: f64, nu: f64, reduced: bool) -> Result<(Vec<f64>, f64)> {
    if xyz.len() < 4 {
        return err("CAX4 braucht 4 Knoten.");
    }
    let mut xy = [[0.0; 2]; 4];
    for i in 0..4 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let d = d_axisym(e, nu)?;
    let n = 8usize;
    let mut ke = vec![0.0; n * n];
    let mut vol = 0.0;
    for (xi, eta, w0) in cax4_gauss(reduced) {
        let (dndx, det, nshp) = quad4_dndx(&xy, xi, eta)?;
        if det <= 0.0 {
            return err("CAX4: negative Jakobideterminante.");
        }
        let r = radius(&xy, &nshp, 4);
        let w = two_pi() * r.max(1e-16) * det * w0;
        let mut b = vec![0.0; 4 * n];
        fill_b_cax(&mut b, 4, &dndx, &nshp, r);
        gemm_bt_d_b(&mut ke, n, &b, 4, &d, w);
        vol += w;
    }
    Ok((ke, vol))
}

pub fn cax8_stiffness(xyz: &[[f64; 3]], e: f64, nu: f64, reduced: bool) -> Result<(Vec<f64>, f64)> {
    if xyz.len() < 8 {
        return err("CAX8 braucht 8 Knoten.");
    }
    let mut xy = [[0.0; 2]; 8];
    for i in 0..8 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let d = d_axisym(e, nu)?;
    let n = 16usize;
    let mut ke = vec![0.0; n * n];
    let mut vol = 0.0;
    for (xi, eta, w0) in cax8_gauss(reduced) {
        let (dndx, det, nshp) = quad8_dndx(&xy, xi, eta)?;
        if det <= 0.0 {
            return err("CAX8: negative Jakobideterminante.");
        }
        let r = radius(&xy, &nshp, 8);
        let w = two_pi() * r.max(1e-16) * det * w0;
        let mut b = vec![0.0; 4 * n];
        fill_b_cax(&mut b, 8, &dndx, &nshp, r);
        gemm_bt_d_b(&mut ke, n, &b, 4, &d, w);
        vol += w;
    }
    Ok((ke, vol))
}

fn voigt_from_cax(s: &[f64]) -> [f64; 6] {
    // εr, εz, εθ, γrz → sxx=σr, syy=σz, szz=σθ, sxy=τrz
    [s[0], s[1], s[2], s[3], 0.0, 0.0]
}

pub fn cax4_nodal_stress(xyz: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<Vec<[f64; 6]>> {
    let mut xy = [[0.0; 2]; 4];
    for i in 0..4 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let d = d_axisym(e, nu)?;
    let n = 8usize;
    let mut out = vec![[0.0; 6]; 4];
    for a in 0..4 {
        let (dndx, _, nshp) = quad4_dndx(&xy, QUAD_XI[a][0], QUAD_XI[a][1])?;
        let r = radius(&xy, &nshp, 4);
        let mut b = vec![0.0; 4 * n];
        fill_b_cax(&mut b, 4, &dndx, &nshp, r);
        let s = sigma_from_b(&b, 4, n, &d, ue);
        out[a] = voigt_from_cax(&s);
    }
    Ok(out)
}

pub fn cax8_nodal_stress(xyz: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<Vec<[f64; 6]>> {
    let mut xy = [[0.0; 2]; 8];
    for i in 0..8 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let d = d_axisym(e, nu)?;
    let n = 16usize;
    let mut out = vec![[0.0; 6]; 8];
    for a in 0..8 {
        let (dndx, _, nshp) = quad8_dndx(&xy, QUAD8_XI[a][0], QUAD8_XI[a][1])?;
        let r = radius(&xy, &nshp, 8);
        let mut b = vec![0.0; 4 * n];
        fill_b_cax(&mut b, 8, &dndx, &nshp, r);
        let s = sigma_from_b(&b, 4, n, &d, ue);
        out[a] = voigt_from_cax(&s);
    }
    Ok(out)
}

pub fn cax4_body_force(xyz: &[[f64; 3]], br: f64, bz: f64, reduced: bool) -> Result<Vec<f64>> {
    let mut xy = [[0.0; 2]; 4];
    for i in 0..4 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let mut fe = vec![0.0; 8];
    for (xi, eta, w0) in cax4_gauss(reduced) {
        let (_, det, nshp) = quad4_dndx(&xy, xi, eta)?;
        if det <= 0.0 {
            continue;
        }
        let r = radius(&xy, &nshp, 4);
        let w = two_pi() * r.max(1e-16) * det * w0;
        for a in 0..4 {
            fe[2 * a] += nshp[a] * br * w;
            fe[2 * a + 1] += nshp[a] * bz * w;
        }
    }
    Ok(fe)
}

pub fn cax8_body_force(xyz: &[[f64; 3]], br: f64, bz: f64, reduced: bool) -> Result<Vec<f64>> {
    if xyz.len() < 8 {
        return err("CAX8 braucht 8 Knoten.");
    }
    let mut xy = [[0.0; 2]; 8];
    for i in 0..8 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let mut fe = vec![0.0; 16];
    for (xi, eta, w0) in cax8_gauss(reduced) {
        let (_, det, nshp) = quad8_dndx(&xy, xi, eta)?;
        if det <= 0.0 {
            continue;
        }
        let r = radius(&xy, &nshp, 8);
        let w = two_pi() * r.max(1e-16) * det * w0;
        for a in 0..8 {
            fe[2 * a] += nshp[a] * br * w;
            fe[2 * a + 1] += nshp[a] * bz * w;
        }
    }
    Ok(fe)
}

fn centrif_at(omega2: f64, origin: [f64; 3], axis_dir: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let mut axis = axis_dir;
    let al = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if al < 1e-18 {
        axis = [0.0, 0.0, 1.0];
    } else {
        axis[0] /= al;
        axis[1] /= al;
        axis[2] /= al;
    }
    let rv = [p[0] - origin[0], p[1] - origin[1], p[2] - origin[2]];
    let proj = rv[0] * axis[0] + rv[1] * axis[1] + rv[2] * axis[2];
    [
        omega2 * (rv[0] - proj * axis[0]),
        omega2 * (rv[1] - proj * axis[1]),
        omega2 * (rv[2] - proj * axis[2]),
    ]
}

fn z_of(xy: &[[f64; 2]], nshp: &[f64], nn: usize) -> f64 {
    let mut z = 0.0;
    for a in 0..nn {
        z += nshp[a] * xy[a][1];
    }
    z
}

/// Centrifugal load with a = ω² r_⊥ evaluated at each Gauss point.
pub fn cax8_centrif_force(
    xyz: &[[f64; 3]],
    rho: f64,
    omega2: f64,
    origin: [f64; 3],
    axis: [f64; 3],
    reduced: bool,
) -> Result<Vec<f64>> {
    if xyz.len() < 8 {
        return err("CAX8 braucht 8 Knoten.");
    }
    let mut xy = [[0.0; 2]; 8];
    for i in 0..8 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let mut fe = vec![0.0; 16];
    for (xi, eta, w0) in cax8_gauss(reduced) {
        let (_, det, nshp) = quad8_dndx(&xy, xi, eta)?;
        if det <= 0.0 {
            continue;
        }
        let r = radius(&xy, &nshp, 8);
        let p = [r, z_of(&xy, &nshp, 8), 0.0];
        let a = centrif_at(omega2, origin, axis, p);
        let w = two_pi() * r.max(1e-16) * det * w0;
        for i in 0..8 {
            fe[2 * i] += nshp[i] * rho * a[0] * w;
            fe[2 * i + 1] += nshp[i] * rho * a[1] * w;
        }
    }
    Ok(fe)
}

pub fn cax4_centrif_force(
    xyz: &[[f64; 3]],
    rho: f64,
    omega2: f64,
    origin: [f64; 3],
    axis: [f64; 3],
    reduced: bool,
) -> Result<Vec<f64>> {
    let mut xy = [[0.0; 2]; 4];
    for i in 0..4.min(xyz.len()) {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let mut fe = vec![0.0; 8];
    for (xi, eta, w0) in cax4_gauss(reduced) {
        let (_, det, nshp) = quad4_dndx(&xy, xi, eta)?;
        if det <= 0.0 {
            continue;
        }
        let r = radius(&xy, &nshp, 4);
        let p = [r, z_of(&xy, &nshp, 4), 0.0];
        let a = centrif_at(omega2, origin, axis, p);
        let w = two_pi() * r.max(1e-16) * det * w0;
        for i in 0..4 {
            fe[2 * i] += nshp[i] * rho * a[0] * w;
            fe[2 * i + 1] += nshp[i] * rho * a[1] * w;
        }
    }
    Ok(fe)
}

pub fn cax_edge_pressure(xyz: &[[f64; 3]], nn: usize, face: i32, p: f64) -> Result<Vec<f64>> {
    // Faces 1-4: edges of the r-z quad. Traction in the outward in-plane normal, 2π r.
    if !(1..=4).contains(&face) {
        return err(format!("Ungültige CAX-Kante P{face}"));
    }
    let edges4 = [[0, 1], [1, 2], [2, 3], [3, 0]];
    let mut fe = vec![0.0; 2 * nn];
    if nn == 3 {
        let edges = [[0, 1], [1, 2], [2, 0]];
        if !(1..=3).contains(&face) {
            return err(format!("Ungültige CAX3-Kante P{face}"));
        }
        let a = edges[(face - 1) as usize][0];
        let b = edges[(face - 1) as usize][1];
        let pa = xyz[a];
        let pb = xyz[b];
        let dr = pb[0] - pa[0];
        let dz = pb[1] - pa[1];
        let len = (dr * dr + dz * dz).sqrt().max(1e-18);
        let nr = dz / len;
        let nz = -dr / len;
        let rmid = 0.5 * (pa[0] + pb[0]).max(0.0);
        let f = -p * len * two_pi() * rmid.max(1e-16) / 2.0;
        fe[2 * a] += f * nr;
        fe[2 * a + 1] += f * nz;
        fe[2 * b] += f * nr;
        fe[2 * b + 1] += f * nz;
        return Ok(fe);
    }
    if nn == 6 {
        // CAX6: 3-node edges 1-2-4, 2-3-5, 3-1-6
        if !(1..=3).contains(&face) {
            return err(format!("Ungültige CAX6-Kante P{face}"));
        }
        let e6 = [[0, 1, 3], [1, 2, 4], [2, 0, 5]];
        let idx = e6[(face - 1) as usize];
        for k in 0..3 {
            let xi = G3[k];
            let w = W3[k];
            let n1 = -0.5 * xi * (1.0 - xi);
            let n2 = 0.5 * xi * (1.0 + xi);
            let n3 = 1.0 - xi * xi;
            let dn1 = xi - 0.5;
            let dn2 = xi + 0.5;
            let dn3 = -2.0 * xi;
            let nshp = [n1, n2, n3];
            let dn = [dn1, dn2, dn3];
            let mut r = 0.0;
            let mut dr = 0.0;
            let mut dz = 0.0;
            for a in 0..3 {
                r += nshp[a] * xyz[idx[a]][0];
                dr += dn[a] * xyz[idx[a]][0];
                dz += dn[a] * xyz[idx[a]][1];
            }
            let jac = (dr * dr + dz * dz).sqrt();
            let nr = dz / jac.max(1e-18);
            let nz = -dr / jac.max(1e-18);
            let coef = -p * jac * two_pi() * r.max(1e-16) * w;
            for a in 0..3 {
                fe[2 * idx[a]] += nshp[a] * coef * nr;
                fe[2 * idx[a] + 1] += nshp[a] * coef * nz;
            }
        }
        return Ok(fe);
    }
    if nn == 4 {
        let a = edges4[(face - 1) as usize][0];
        let b = edges4[(face - 1) as usize][1];
        let pa = xyz[a];
        let pb = xyz[b];
        let dr = pb[0] - pa[0];
        let dz = pb[1] - pa[1];
        let len = (dr * dr + dz * dz).sqrt().max(1e-18);
        // outward = rotate edge 90° (assuming CCW): (dz, -dr)
        let nr = dz / len;
        let nz = -dr / len;
        let rmid = 0.5 * (pa[0] + pb[0]).max(0.0);
        let f = -p * len * two_pi() * rmid.max(1e-16) / 2.0;
        fe[2 * a] += f * nr;
        fe[2 * a + 1] += f * nz;
        fe[2 * b] += f * nr;
        fe[2 * b + 1] += f * nz;
    } else {
        // 8-node: 3-node edge, Gauss on [-1,1]
        // face 1: 1-2-5, face 2: 2-3-6, face 3: 3-4-7, face 4: 4-1-8
        let e8 = [[0, 1, 4], [1, 2, 5], [2, 3, 6], [3, 0, 7]];
        let idx = e8[(face - 1) as usize];
        for k in 0..3 {
            let xi = G3[k];
            let w = W3[k];
            let n1 = -0.5 * xi * (1.0 - xi);
            let n2 = 0.5 * xi * (1.0 + xi);
            let n3 = 1.0 - xi * xi;
            let dn1 = xi - 0.5;
            let dn2 = xi + 0.5;
            let dn3 = -2.0 * xi;
            let nshp = [n1, n2, n3];
            let dn = [dn1, dn2, dn3];
            let mut r = 0.0;
            let mut z = 0.0;
            let mut dr = 0.0;
            let mut dz = 0.0;
            for a in 0..3 {
                r += nshp[a] * xyz[idx[a]][0];
                z += nshp[a] * xyz[idx[a]][1];
                dr += dn[a] * xyz[idx[a]][0];
                dz += dn[a] * xyz[idx[a]][1];
            }
            let _ = z;
            let jac = (dr * dr + dz * dz).sqrt();
            let nr = dz / jac.max(1e-18);
            let nz = -dr / jac.max(1e-18);
            let coef = -p * jac * two_pi() * r.max(1e-16) * w;
            for a in 0..3 {
                fe[2 * idx[a]] += nshp[a] * coef * nr;
                fe[2 * idx[a] + 1] += nshp[a] * coef * nz;
            }
        }
    }
    Ok(fe)
}

fn tri_gauss3() -> [(f64, f64, f64); 3] {
    [
        (1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0),
        (2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0),
        (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0),
    ]
}

fn tri3_dndx(xy: &[[f64; 2]]) -> Result<([[f64; 2]; 3], f64)> {
    let two_a = xy[0][0] * (xy[1][1] - xy[2][1])
        + xy[1][0] * (xy[2][1] - xy[0][1])
        + xy[2][0] * (xy[0][1] - xy[1][1]);
    if two_a <= 0.0 {
        return err("CAX3: nicht-positive Fläche (Knotenreihenfolge).");
    }
    let mut dndx = [[0.0; 2]; 3];
    dndx[0] = [(xy[1][1] - xy[2][1]) / two_a, (xy[2][0] - xy[1][0]) / two_a];
    dndx[1] = [(xy[2][1] - xy[0][1]) / two_a, (xy[0][0] - xy[2][0]) / two_a];
    dndx[2] = [(xy[0][1] - xy[1][1]) / two_a, (xy[1][0] - xy[0][0]) / two_a];
    Ok((dndx, two_a))
}

pub fn cax3_stiffness(xyz: &[[f64; 3]], e: f64, nu: f64) -> Result<(Vec<f64>, f64)> {
    if xyz.len() < 3 {
        return err("CAX3 braucht 3 Knoten.");
    }
    let mut xy = [[0.0; 2]; 3];
    for i in 0..3 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let d = d_axisym(e, nu)?;
    let n = 6usize;
    let mut ke = vec![0.0; n * n];
    let mut vol = 0.0;
    let (dndx, two_a) = tri3_dndx(&xy)?;
    for (xi, eta, w0) in tri_gauss3() {
        let nshp = [1.0 - xi - eta, xi, eta];
        let r = radius(&xy, &nshp, 3);
        let w = two_pi() * r.max(1e-16) * two_a * w0;
        let mut b = vec![0.0; 4 * n];
        fill_b_cax(&mut b, 3, &dndx, &nshp, r);
        gemm_bt_d_b(&mut ke, n, &b, 4, &d, w);
        vol += w;
    }
    Ok((ke, vol))
}

pub fn cax6_stiffness(xyz: &[[f64; 3]], e: f64, nu: f64) -> Result<(Vec<f64>, f64)> {
    if xyz.len() < 6 {
        return err("CAX6 braucht 6 Knoten.");
    }
    let mut xy = [[0.0; 2]; 6];
    for i in 0..6 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let d = d_axisym(e, nu)?;
    let n = 12usize;
    let mut ke = vec![0.0; n * n];
    let mut vol = 0.0;
    for (xi, eta, w0) in tri_gauss3() {
        let (dndx, det, nshp) = crate::quadratic::tri6_dndx(&xy, xi, eta)?;
        if det <= 0.0 {
            return err("CAX6: negative Jakobideterminante.");
        }
        let r = radius(&xy, &nshp, 6);
        let w = two_pi() * r.max(1e-16) * det * w0;
        let mut b = vec![0.0; 4 * n];
        fill_b_cax(&mut b, 6, &dndx, &nshp, r);
        gemm_bt_d_b(&mut ke, n, &b, 4, &d, w);
        vol += w;
    }
    Ok((ke, vol))
}

pub fn cax3_nodal_stress(xyz: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<Vec<[f64; 6]>> {
    let mut xy = [[0.0; 2]; 3];
    for i in 0..3 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let d = d_axisym(e, nu)?;
    let n = 6usize;
    let (dndx, _) = tri3_dndx(&xy)?;
    let mut out = vec![[0.0; 6]; 3];
    for a in 0..3 {
        let mut nshp = [0.0; 3];
        nshp[a] = 1.0;
        let r = radius(&xy, &nshp, 3);
        let mut b = vec![0.0; 4 * n];
        fill_b_cax(&mut b, 3, &dndx, &nshp, r);
        let s = sigma_from_b(&b, 4, n, &d, ue);
        out[a] = voigt_from_cax(&s);
    }
    Ok(out)
}

pub fn cax6_nodal_stress(xyz: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<Vec<[f64; 6]>> {
    let mut xy = [[0.0; 2]; 6];
    for i in 0..6 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let d = d_axisym(e, nu)?;
    let n = 12usize;
    let corners = [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.5, 0.0],
        [0.5, 0.5],
        [0.0, 0.5],
    ];
    let mut out = vec![[0.0; 6]; 6];
    for a in 0..6 {
        let (dndx, _, nshp) = crate::quadratic::tri6_dndx(&xy, corners[a][0], corners[a][1])?;
        let r = radius(&xy, &nshp, 6);
        let mut b = vec![0.0; 4 * n];
        fill_b_cax(&mut b, 6, &dndx, &nshp, r);
        let s = sigma_from_b(&b, 4, n, &d, ue);
        out[a] = voigt_from_cax(&s);
    }
    Ok(out)
}

pub fn cax3_body_force(xyz: &[[f64; 3]], br: f64, bz: f64) -> Result<Vec<f64>> {
    let mut xy = [[0.0; 2]; 3];
    for i in 0..3 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let (_, two_a) = tri3_dndx(&xy)?;
    let mut fe = vec![0.0; 6];
    for (xi, eta, w0) in tri_gauss3() {
        let nshp = [1.0 - xi - eta, xi, eta];
        let r = radius(&xy, &nshp, 3);
        let w = two_pi() * r.max(1e-16) * two_a * w0;
        for a in 0..3 {
            fe[2 * a] += nshp[a] * br * w;
            fe[2 * a + 1] += nshp[a] * bz * w;
        }
    }
    Ok(fe)
}

pub fn cax6_body_force(xyz: &[[f64; 3]], br: f64, bz: f64) -> Result<Vec<f64>> {
    let mut xy = [[0.0; 2]; 6];
    for i in 0..6 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let mut fe = vec![0.0; 12];
    for (xi, eta, w0) in tri_gauss3() {
        let (_, det, nshp) = crate::quadratic::tri6_dndx(&xy, xi, eta)?;
        let r = radius(&xy, &nshp, 6);
        let w = two_pi() * r.max(1e-16) * det * w0;
        for a in 0..6 {
            fe[2 * a] += nshp[a] * br * w;
            fe[2 * a + 1] += nshp[a] * bz * w;
        }
    }
    Ok(fe)
}
