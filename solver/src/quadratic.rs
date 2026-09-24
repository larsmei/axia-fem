//! Quadratic continuum elements: C3D20(R), C3D10/C3D10T, CPS8/CPE8(R), CPS6/CPE6.
//! Node numbering matches Abaqus / CalculiX.

use crate::elem::{
    d_iso_3d, d_plane_strain, d_plane_stress, fill_b2, fill_b3, gemm_bt_d_b, invert2, invert3,
    sigma_from_b, G2,
};
use crate::error::{err, Result};

pub(crate) const G3: [f64; 3] = [-0.7745966692414834, 0.0, 0.7745966692414834];
pub(crate) const W3: [f64; 3] = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];

pub(crate) const HEX20_XI: [[f64; 3]; 20] = [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
    [0.0, -1.0, -1.0],
    [1.0, 0.0, -1.0],
    [0.0, 1.0, -1.0],
    [-1.0, 0.0, -1.0],
    [0.0, -1.0, 1.0],
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 1.0],
    [-1.0, 0.0, 1.0],
    [-1.0, -1.0, 0.0],
    [1.0, -1.0, 0.0],
    [1.0, 1.0, 0.0],
    [-1.0, 1.0, 0.0],
];

fn hex20_shape(xi: f64, eta: f64, zeta: f64) -> ([f64; 20], [[f64; 3]; 20]) {
    let mut n = [0.0; 20];
    let mut dn = [[0.0; 3]; 20];
    for i in 0..8 {
        let x = HEX20_XI[i][0];
        let e = HEX20_XI[i][1];
        let z = HEX20_XI[i][2];
        let a = 1.0 + x * xi;
        let b = 1.0 + e * eta;
        let c = 1.0 + z * zeta;
        n[i] = 0.125 * a * b * c * (x * xi + e * eta + z * zeta - 2.0);
        dn[i][0] = 0.125 * x * b * c * (2.0 * x * xi + e * eta + z * zeta - 1.0);
        dn[i][1] = 0.125 * e * a * c * (x * xi + 2.0 * e * eta + z * zeta - 1.0);
        dn[i][2] = 0.125 * z * a * b * (x * xi + e * eta + 2.0 * z * zeta - 1.0);
    }
    for i in 8..20 {
        let x = HEX20_XI[i][0];
        let e = HEX20_XI[i][1];
        let z = HEX20_XI[i][2];
        if x.abs() < 0.5 {
            n[i] = 0.25 * (1.0 - xi * xi) * (1.0 + e * eta) * (1.0 + z * zeta);
            dn[i][0] = 0.25 * (-2.0 * xi) * (1.0 + e * eta) * (1.0 + z * zeta);
            dn[i][1] = 0.25 * (1.0 - xi * xi) * e * (1.0 + z * zeta);
            dn[i][2] = 0.25 * (1.0 - xi * xi) * (1.0 + e * eta) * z;
        } else if e.abs() < 0.5 {
            n[i] = 0.25 * (1.0 - eta * eta) * (1.0 + x * xi) * (1.0 + z * zeta);
            dn[i][0] = 0.25 * (1.0 - eta * eta) * x * (1.0 + z * zeta);
            dn[i][1] = 0.25 * (-2.0 * eta) * (1.0 + x * xi) * (1.0 + z * zeta);
            dn[i][2] = 0.25 * (1.0 - eta * eta) * (1.0 + x * xi) * z;
        } else {
            n[i] = 0.25 * (1.0 - zeta * zeta) * (1.0 + x * xi) * (1.0 + e * eta);
            dn[i][0] = 0.25 * (1.0 - zeta * zeta) * x * (1.0 + e * eta);
            dn[i][1] = 0.25 * (1.0 - zeta * zeta) * (1.0 + x * xi) * e;
            dn[i][2] = 0.25 * (-2.0 * zeta) * (1.0 + x * xi) * (1.0 + e * eta);
        }
    }
    (n, dn)
}

pub(crate) fn hex20_dndx(
    xyz: &[[f64; 3]],
    xi: f64,
    eta: f64,
    zeta: f64,
) -> Result<([[f64; 3]; 20], f64, [f64; 20])> {
    let (n, dn) = hex20_shape(xi, eta, zeta);
    let mut j = [[0.0; 3]; 3];
    for a in 0..20 {
        for p in 0..3 {
            for q in 0..3 {
                j[q][p] += dn[a][p] * xyz[a][q];
            }
        }
    }
    let (inv, det) = invert3(j)?;
    let mut dndx = [[0.0; 3]; 20];
    for a in 0..20 {
        for i in 0..3 {
            dndx[a][i] = inv[0][i] * dn[a][0] + inv[1][i] * dn[a][1] + inv[2][i] * dn[a][2];
        }
    }
    Ok((dndx, det, n))
}

pub(crate) fn hex_gauss(reduced: bool) -> Vec<(f64, f64, f64, f64)> {
    let mut o = Vec::new();
    if reduced {
        for &xi in &[-G2, G2] {
            for &eta in &[-G2, G2] {
                for &zeta in &[-G2, G2] {
                    o.push((xi, eta, zeta, 1.0));
                }
            }
        }
    } else {
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    o.push((G3[i], G3[j], G3[k], W3[i] * W3[j] * W3[k]));
                }
            }
        }
    }
    o
}

pub fn hex20_stiffness(
    xyz: &[[f64; 3]],
    e: f64,
    nu: f64,
    reduced: bool,
) -> Result<(Vec<f64>, f64)> {
    if xyz.len() < 20 {
        return err("C3D20 braucht 20 Knoten.");
    }
    let d = d_iso_3d(e, nu)?;
    let n = 60usize;
    let mut ke = vec![0.0; n * n];
    let mut vol = 0.0;
    for (xi, eta, zeta, w) in hex_gauss(reduced) {
        let (dndx, det, _) = hex20_dndx(xyz, xi, eta, zeta)?;
        if det <= 0.0 {
            return err("C3D20: negative Jakobideterminante.");
        }
        let mut b = vec![0.0; 6 * n];
        fill_b3(&mut b, 20, &dndx);
        gemm_bt_d_b(&mut ke, n, &b, 6, &d, w * det);
        vol += w * det;
    }
    Ok((ke, vol))
}

pub fn hex20_nodal_stress(
    xyz: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    reduced: bool,
) -> Result<Vec<[f64; 6]>> {
    // Gauss-point stresses, then Lagrange extrapolation to the 20 nodes.
    // Evaluating B at ξ=±1 on a serendipity hex wildly overshoots peaks
    // (bolted-joint vmMax was ~2.4× CalculiX).
    let d = d_iso_3d(e, nu)?;
    let n = 60usize;
    let gps = hex_gauss(reduced);
    let mut gsig = Vec::with_capacity(gps.len());
    for &(xi, eta, zeta, _) in &gps {
        let (dndx, _, _) = hex20_dndx(xyz, xi, eta, zeta)?;
        let mut b = vec![0.0; 6 * n];
        fill_b3(&mut b, 20, &dndx);
        let s = sigma_from_b(&b, 6, n, &d, ue);
        let mut six = [0.0; 6];
        six.copy_from_slice(&s);
        gsig.push(six);
    }
    Ok(hex20_extrapolate(&gsig, reduced))
}

pub(crate) fn hex20_extrapolate(gsig: &[[f64; 6]], reduced: bool) -> Vec<[f64; 6]> {
    let pts: Vec<f64> = if reduced { vec![-G2, G2] } else { G3.to_vec() };
    let n1 = pts.len();
    let mut out = vec![[0.0; 6]; 20];
    for a in 0..20 {
        let lx = lagrange1d(&pts, HEX20_XI[a][0]);
        let ly = lagrange1d(&pts, HEX20_XI[a][1]);
        let lz = lagrange1d(&pts, HEX20_XI[a][2]);
        let mut idx = 0usize;
        for i in 0..n1 {
            for j in 0..n1 {
                for k in 0..n1 {
                    let w = lx[i] * ly[j] * lz[k];
                    for c in 0..6 {
                        out[a][c] += w * gsig[idx][c];
                    }
                    idx += 1;
                }
            }
        }
    }
    out
}

fn lagrange1d(pts: &[f64], x: f64) -> Vec<f64> {
    let n = pts.len();
    let mut w = vec![1.0; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            w[i] *= (x - pts[j]) / (pts[i] - pts[j]);
        }
    }
    w
}

pub fn hex20_body_force(
    xyz: &[[f64; 3]],
    bx: f64,
    by: f64,
    bz: f64,
    reduced: bool,
) -> Result<Vec<f64>> {
    let mut fe = vec![0.0; 60];
    for (xi, eta, zeta, w) in hex_gauss(reduced) {
        let (_, det, n) = hex20_dndx(xyz, xi, eta, zeta)?;
        for a in 0..20 {
            fe[3 * a] += n[a] * bx * w * det;
            fe[3 * a + 1] += n[a] * by * w * det;
            fe[3 * a + 2] += n[a] * bz * w * det;
        }
    }
    Ok(fe)
}

/// CalculiX P1..P6 on C3D20. 8-node faces: 4 corners + 4 midsides.
pub(crate) const HEX20_FACE: [[usize; 8]; 6] = [
    [0, 1, 2, 3, 8, 9, 10, 11],
    [4, 7, 6, 5, 15, 14, 13, 12],
    [0, 4, 5, 1, 16, 12, 17, 8],
    [1, 5, 6, 2, 17, 13, 18, 9],
    [2, 6, 7, 3, 18, 14, 19, 10],
    [3, 7, 4, 0, 19, 15, 16, 11],
];

pub fn hex20_face_pressure(xyz: &[[f64; 3]], face: i32, p: f64) -> Result<Vec<f64>> {
    let faces = HEX20_FACE;
    if !(1..=6).contains(&face) {
        return err(format!("Ungültige C3D20-Fläche P{face}"));
    }
    let fi = (face - 1) as usize;
    let mut fe = vec![0.0; 60];
    let mut gps = Vec::new();
    for i in 0..3 {
        for j in 0..3 {
            gps.push((G3[i], G3[j], W3[i] * W3[j]));
        }
    }
    for (xi, eta, w) in gps {
        let (nshp, dn) = quad8_shape(xi, eta);
        let mut rxi = [0.0; 3];
        let mut reta = [0.0; 3];
        for a in 0..8 {
            let q = xyz[faces[fi][a]];
            for k in 0..3 {
                rxi[k] += dn[a][0] * q[k];
                reta[k] += dn[a][1] * q[k];
            }
        }
        let nx = rxi[1] * reta[2] - rxi[2] * reta[1];
        let ny = rxi[2] * reta[0] - rxi[0] * reta[2];
        let nz = rxi[0] * reta[1] - rxi[1] * reta[0];
        let tx = -p * nx * w;
        let ty = -p * ny * w;
        let tz = -p * nz * w;
        for a in 0..8 {
            let gd = 3 * faces[fi][a];
            fe[gd] += nshp[a] * tx;
            fe[gd + 1] += nshp[a] * ty;
            fe[gd + 2] += nshp[a] * tz;
        }
    }
    Ok(fe)
}

pub(crate) const QUAD8_XI: [[f64; 2]; 8] = [
    [-1.0, -1.0],
    [1.0, -1.0],
    [1.0, 1.0],
    [-1.0, 1.0],
    [0.0, -1.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [-1.0, 0.0],
];

pub(crate) fn quad8_shape(xi: f64, eta: f64) -> ([f64; 8], [[f64; 2]; 8]) {
    let mut n = [0.0; 8];
    let mut dn = [[0.0; 2]; 8];
    for i in 0..4 {
        let x = QUAD8_XI[i][0];
        let e = QUAD8_XI[i][1];
        let a = 1.0 + x * xi;
        let b = 1.0 + e * eta;
        n[i] = 0.25 * a * b * (x * xi + e * eta - 1.0);
        dn[i][0] = 0.25 * x * b * (2.0 * x * xi + e * eta);
        dn[i][1] = 0.25 * e * a * (x * xi + 2.0 * e * eta);
    }
    for i in 4..8 {
        let x = QUAD8_XI[i][0];
        let e = QUAD8_XI[i][1];
        if x.abs() < 0.5 {
            n[i] = 0.5 * (1.0 - xi * xi) * (1.0 + e * eta);
            dn[i][0] = -xi * (1.0 + e * eta);
            dn[i][1] = 0.5 * (1.0 - xi * xi) * e;
        } else {
            n[i] = 0.5 * (1.0 - eta * eta) * (1.0 + x * xi);
            dn[i][0] = 0.5 * (1.0 - eta * eta) * x;
            dn[i][1] = -eta * (1.0 + x * xi);
        }
    }
    (n, dn)
}

pub(crate) fn quad8_dndx(
    xy: &[[f64; 2]],
    xi: f64,
    eta: f64,
) -> Result<([[f64; 2]; 8], f64, [f64; 8])> {
    let (n, dn) = quad8_shape(xi, eta);
    let mut j = [[0.0; 2]; 2];
    for a in 0..8 {
        j[0][0] += dn[a][0] * xy[a][0];
        j[0][1] += dn[a][1] * xy[a][0];
        j[1][0] += dn[a][0] * xy[a][1];
        j[1][1] += dn[a][1] * xy[a][1];
    }
    let (inv, det) = invert2(j)?;
    let mut dndx = [[0.0; 2]; 8];
    for a in 0..8 {
        dndx[a][0] = inv[0][0] * dn[a][0] + inv[1][0] * dn[a][1];
        dndx[a][1] = inv[0][1] * dn[a][0] + inv[1][1] * dn[a][1];
    }
    Ok((dndx, det, n))
}

fn quad_gauss(reduced: bool) -> Vec<(f64, f64, f64)> {
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

pub fn quad8_stiffness(
    xy: &[[f64; 2]],
    e: f64,
    nu: f64,
    t: f64,
    plane_strain: bool,
    reduced: bool,
) -> Result<(Vec<f64>, f64)> {
    let d = if plane_strain {
        d_plane_strain(e, nu)?
    } else {
        d_plane_stress(e, nu)?
    };
    let n = 16usize;
    let mut ke = vec![0.0; n * n];
    let mut area = 0.0;
    for (xi, eta, w) in quad_gauss(reduced) {
        let (dndx, det, _) = quad8_dndx(xy, xi, eta)?;
        if det <= 0.0 {
            return err("CPS8/CPE8: negative Jakobideterminante.");
        }
        let mut b = vec![0.0; 3 * n];
        fill_b2(&mut b, 8, &dndx);
        gemm_bt_d_b(&mut ke, n, &b, 3, &d, t * w * det);
        area += w * det;
    }
    Ok((ke, area))
}

pub fn quad8_nodal_stress(
    xyz: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    plane_strain: bool,
) -> Result<Vec<[f64; 6]>> {
    let mut xy = [[0.0; 2]; 8];
    for i in 0..8 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let d = if plane_strain {
        d_plane_strain(e, nu)?
    } else {
        d_plane_stress(e, nu)?
    };
    let n = 16usize;
    let mut out = vec![[0.0; 6]; 8];
    for a in 0..8 {
        let (dndx, _, _) = quad8_dndx(&xy, QUAD8_XI[a][0], QUAD8_XI[a][1])?;
        let mut b = vec![0.0; 3 * n];
        fill_b2(&mut b, 8, &dndx);
        let s = sigma_from_b(&b, 3, n, &d, ue);
        out[a][0] = s[0];
        out[a][1] = s[1];
        out[a][3] = s[2];
        if plane_strain {
            // σz = ν(σx+σy)
            // recovered from plane strain D; s only has xx,yy,xy
        }
    }
    Ok(out)
}

pub fn quad8_body_force(
    xy: &[[f64; 2]],
    bx: f64,
    by: f64,
    t: f64,
    reduced: bool,
) -> Result<Vec<f64>> {
    let mut fe = vec![0.0; 16];
    for (xi, eta, w) in quad_gauss(reduced) {
        let (_, det, n) = quad8_dndx(xy, xi, eta)?;
        for a in 0..8 {
            fe[2 * a] += n[a] * bx * w * det * t;
            fe[2 * a + 1] += n[a] * by * w * det * t;
        }
    }
    Ok(fe)
}

pub fn quad8_edge_pressure(xy: &[[f64; 2]], face: i32, p: f64, thickness: f64) -> Result<Vec<f64>> {
    // P1..P4: edge with two corners + midside. Indices into 8-node connectivity.
    let edges: [[usize; 3]; 4] = [[0, 1, 4], [1, 2, 5], [2, 3, 6], [3, 0, 7]];
    if !(1..=4).contains(&face) {
        return err(format!("Ungültige CPS8-Kante P{face}"));
    }
    let e = (face - 1) as usize;
    let mut fe = vec![0.0; 16];
    // 3-node quadratic edge, 2-point Gauss in ξ∈[-1,1]
    for &(xi, w) in &[(-G2, 1.0), (G2, 1.0)] {
        let n1 = 0.5 * xi * (xi - 1.0);
        let n2 = 0.5 * xi * (xi + 1.0);
        let n3 = 1.0 - xi * xi;
        let dn1 = xi - 0.5;
        let dn2 = xi + 0.5;
        let dn3 = -2.0 * xi;
        let a = edges[e][0];
        let b = edges[e][1];
        let m = edges[e][2];
        let dx = dn1 * xy[a][0] + dn2 * xy[b][0] + dn3 * xy[m][0];
        let dy = dn1 * xy[a][1] + dn2 * xy[b][1] + dn3 * xy[m][1];
        // inward for CCW: (-dy, dx)
        let fx = p * thickness * (-dy) * w;
        let fy = p * thickness * dx * w;
        fe[2 * a] += n1 * fx;
        fe[2 * a + 1] += n1 * fy;
        fe[2 * b] += n2 * fx;
        fe[2 * b + 1] += n2 * fy;
        fe[2 * m] += n3 * fx;
        fe[2 * m + 1] += n3 * fy;
    }
    Ok(fe)
}

/// CPS6/CPE6 edge pressure. P1 = edge 1-2, P2 = 2-3, P3 = 3-1.
/// Positive pressure points into a CCW element.
pub fn tri6_edge_pressure(xy: &[[f64; 2]], face: i32, p: f64, thickness: f64) -> Result<Vec<f64>> {
    let edges: [[usize; 3]; 3] = [[0, 1, 3], [1, 2, 4], [2, 0, 5]];
    if !(1..=3).contains(&face) {
        return err(format!("Ungültige CPS6-Kante P{face}"));
    }
    let e = (face - 1) as usize;
    let mut fe = vec![0.0; 12];
    for &(xi, w) in &[(-G2, 1.0), (G2, 1.0)] {
        let n1 = 0.5 * xi * (xi - 1.0);
        let n2 = 0.5 * xi * (xi + 1.0);
        let n3 = 1.0 - xi * xi;
        let dn1 = xi - 0.5;
        let dn2 = xi + 0.5;
        let dn3 = -2.0 * xi;
        let a = edges[e][0];
        let b = edges[e][1];
        let m = edges[e][2];
        let dx = dn1 * xy[a][0] + dn2 * xy[b][0] + dn3 * xy[m][0];
        let dy = dn1 * xy[a][1] + dn2 * xy[b][1] + dn3 * xy[m][1];
        let fx = p * thickness * (-dy) * w;
        let fy = p * thickness * dx * w;
        fe[2 * a] += n1 * fx;
        fe[2 * a + 1] += n1 * fy;
        fe[2 * b] += n2 * fx;
        fe[2 * b + 1] += n2 * fy;
        fe[2 * m] += n3 * fx;
        fe[2 * m + 1] += n3 * fy;
    }
    Ok(fe)
}

// ---- C3D10 / CPS6 ----------------------------------------------------------

fn tet10_shape(r: f64, s: f64, t: f64) -> ([f64; 10], [[f64; 3]; 10]) {
    let l1 = 1.0 - r - s - t;
    let l2 = r;
    let l3 = s;
    let l4 = t;
    let n = [
        l1 * (2.0 * l1 - 1.0),
        l2 * (2.0 * l2 - 1.0),
        l3 * (2.0 * l3 - 1.0),
        l4 * (2.0 * l4 - 1.0),
        4.0 * l1 * l2,
        4.0 * l2 * l3,
        4.0 * l3 * l1,
        4.0 * l1 * l4,
        4.0 * l2 * l4,
        4.0 * l3 * l4,
    ];
    let d1 = 4.0 * l1 - 1.0;
    let d2 = 4.0 * l2 - 1.0;
    let d3 = 4.0 * l3 - 1.0;
    let d4 = 4.0 * l4 - 1.0;
    let dn = [
        [-d1, -d1, -d1],
        [d2, 0.0, 0.0],
        [0.0, d3, 0.0],
        [0.0, 0.0, d4],
        [4.0 * (l1 - l2), -4.0 * l2, -4.0 * l2],
        [4.0 * l3, 4.0 * l2, 0.0],
        [-4.0 * l3, 4.0 * (l1 - l3), -4.0 * l3],
        [-4.0 * l4, -4.0 * l4, 4.0 * (l1 - l4)],
        [4.0 * l4, 0.0, 4.0 * l2],
        [0.0, 4.0 * l4, 4.0 * l3],
    ];
    (n, dn)
}

pub(crate) fn tet10_dndx(
    xyz: &[[f64; 3]],
    r: f64,
    s: f64,
    t: f64,
) -> Result<([[f64; 3]; 10], f64, [f64; 10])> {
    let (n, dn) = tet10_shape(r, s, t);
    let mut j = [[0.0; 3]; 3];
    for a in 0..10 {
        for p in 0..3 {
            for q in 0..3 {
                j[q][p] += dn[a][p] * xyz[a][q];
            }
        }
    }
    let (inv, det) = invert3(j)?;
    let mut dndx = [[0.0; 3]; 10];
    for a in 0..10 {
        for i in 0..3 {
            dndx[a][i] = inv[0][i] * dn[a][0] + inv[1][i] * dn[a][1] + inv[2][i] * dn[a][2];
        }
    }
    Ok((dndx, det, n))
}

pub fn tet10_stiffness(xyz: &[[f64; 3]], e: f64, nu: f64) -> Result<(Vec<f64>, f64)> {
    let d = d_iso_3d(e, nu)?;
    let n = 30usize;
    let mut ke = vec![0.0; n * n];
    let mut vol = 0.0;
    let a = 0.5854101966249685;
    let b = 0.1381966011250105;
    let w = 1.0 / 24.0;
    let pts = [[b, b, b], [a, b, b], [b, a, b], [b, b, a]];
    for p in &pts {
        let (dndx, det, _) = tet10_dndx(xyz, p[0], p[1], p[2])?;
        if det <= 0.0 {
            return err("C3D10: negative Jakobideterminante (Knotenreihenfolge).");
        }
        let mut bb = vec![0.0; 6 * n];
        fill_b3(&mut bb, 10, &dndx);
        gemm_bt_d_b(&mut ke, n, &bb, 6, &d, w * det);
        vol += w * det;
    }
    Ok((ke, vol))
}

fn fill_b_vol(b: &mut [f64], nnode: usize, dndx: &[[f64; 3]]) {
    let n = 3 * nnode;
    for i in 0..nnode {
        let c = 3 * i;
        let dx = dndx[i][0] / 3.0;
        let dy = dndx[i][1] / 3.0;
        let dz = dndx[i][2] / 3.0;
        for row in 0..3 {
            b[row * n + c] += dx;
            b[row * n + c + 1] += dy;
            b[row * n + c + 2] += dz;
        }
    }
}

/// C3D10T: tet10 with mean-dilatation (B-bar) volumetric strain.
pub fn tet10t_stiffness(xyz: &[[f64; 3]], e: f64, nu: f64) -> Result<(Vec<f64>, f64)> {
    let d = d_iso_3d(e, nu)?;
    let n = 30usize;
    let mut ke = vec![0.0; n * n];
    let mut vol = 0.0;
    let a = 0.5854101966249685;
    let b = 0.1381966011250105;
    let w = 1.0 / 24.0;
    let pts = [[b, b, b], [a, b, b], [b, a, b], [b, b, a]];
    let (dndx0, det0, _) = tet10_dndx(xyz, 0.25, 0.25, 0.25)?;
    if det0 <= 0.0 {
        return err("C3D10T: negative Jakobideterminante.");
    }
    let mut bvol0 = vec![0.0; 6 * n];
    fill_b_vol(&mut bvol0, 10, &dndx0);
    for p in &pts {
        let (dndx, det, _) = tet10_dndx(xyz, p[0], p[1], p[2])?;
        if det <= 0.0 {
            return err("C3D10T: negative Jakobideterminante (Knotenreihenfolge).");
        }
        let mut bb = vec![0.0; 6 * n];
        fill_b3(&mut bb, 10, &dndx);
        let mut bvol = vec![0.0; 6 * n];
        fill_b_vol(&mut bvol, 10, &dndx);
        for i in 0..bb.len() {
            bb[i] += bvol0[i] - bvol[i];
        }
        gemm_bt_d_b(&mut ke, n, &bb, 6, &d, w * det);
        vol += w * det;
    }
    Ok((ke, vol))
}

pub fn tet10_nodal_stress(xyz: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<Vec<[f64; 6]>> {
    let d = d_iso_3d(e, nu)?;
    let n = 30usize;
    let rst = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    let mut out = vec![[0.0; 6]; 10];
    for a in 0..10 {
        let (dndx, _, _) = tet10_dndx(xyz, rst[a][0], rst[a][1], rst[a][2])?;
        let mut b = vec![0.0; 6 * n];
        fill_b3(&mut b, 10, &dndx);
        let s = sigma_from_b(&b, 6, n, &d, ue);
        out[a].copy_from_slice(&s);
    }
    Ok(out)
}

pub fn tet10_body_force(xyz: &[[f64; 3]], bx: f64, by: f64, bz: f64) -> Result<Vec<f64>> {
    let mut fe = vec![0.0; 30];
    let a = 0.5854101966249685;
    let b = 0.1381966011250105;
    let w = 1.0 / 24.0;
    let pts = [[b, b, b], [a, b, b], [b, a, b], [b, b, a]];
    for p in &pts {
        let (_, det, n) = tet10_dndx(xyz, p[0], p[1], p[2])?;
        for i in 0..10 {
            fe[3 * i] += n[i] * bx * w * det;
            fe[3 * i + 1] += n[i] * by * w * det;
            fe[3 * i + 2] += n[i] * bz * w * det;
        }
    }
    Ok(fe)
}

pub(crate) fn tri6_shape(xi: f64, eta: f64) -> ([f64; 6], [[f64; 2]; 6]) {
    let l1 = 1.0 - xi - eta;
    let l2 = xi;
    let l3 = eta;
    let n = [
        l1 * (2.0 * l1 - 1.0),
        l2 * (2.0 * l2 - 1.0),
        l3 * (2.0 * l3 - 1.0),
        4.0 * l1 * l2,
        4.0 * l2 * l3,
        4.0 * l3 * l1,
    ];
    let d1 = 4.0 * l1 - 1.0;
    let d2 = 4.0 * l2 - 1.0;
    let d3 = 4.0 * l3 - 1.0;
    let mut dn = [[0.0; 2]; 6];
    dn[0] = [-d1, -d1];
    dn[1] = [d2, 0.0];
    dn[2] = [0.0, d3];
    dn[3] = [4.0 * (l1 - l2), -4.0 * l2];
    dn[4] = [4.0 * l3, 4.0 * l2];
    dn[5] = [-4.0 * l3, 4.0 * (l1 - l3)];
    (n, dn)
}

pub(crate) fn tri6_dndx(
    xy: &[[f64; 2]],
    xi: f64,
    eta: f64,
) -> Result<([[f64; 2]; 6], f64, [f64; 6])> {
    let (n, dn) = tri6_shape(xi, eta);
    let mut j = [[0.0; 2]; 2];
    for a in 0..6 {
        j[0][0] += dn[a][0] * xy[a][0];
        j[0][1] += dn[a][1] * xy[a][0];
        j[1][0] += dn[a][0] * xy[a][1];
        j[1][1] += dn[a][1] * xy[a][1];
    }
    let (inv, det) = invert2(j)?;
    let mut dndx = [[0.0; 2]; 6];
    for a in 0..6 {
        dndx[a][0] = inv[0][0] * dn[a][0] + inv[1][0] * dn[a][1];
        dndx[a][1] = inv[0][1] * dn[a][0] + inv[1][1] * dn[a][1];
    }
    Ok((dndx, det, n))
}

pub fn tri6_stiffness(
    xy: &[[f64; 2]],
    e: f64,
    nu: f64,
    t: f64,
    plane_strain: bool,
) -> Result<(Vec<f64>, f64)> {
    let d = if plane_strain {
        d_plane_strain(e, nu)?
    } else {
        d_plane_stress(e, nu)?
    };
    let n = 12usize;
    let mut ke = vec![0.0; n * n];
    let mut area = 0.0;
    // 3-point Hammer, parent area 1/2
    let pts = [
        [1.0 / 6.0, 1.0 / 6.0],
        [2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0],
    ];
    let w = 1.0 / 6.0;
    for p in &pts {
        let (dndx, det, _) = tri6_dndx(xy, p[0], p[1])?;
        if det <= 0.0 {
            return err("CPS6/CPE6: negative Jakobideterminante.");
        }
        let mut b = vec![0.0; 3 * n];
        fill_b2(&mut b, 6, &dndx);
        gemm_bt_d_b(&mut ke, n, &b, 3, &d, t * w * det);
        area += w * det;
    }
    Ok((ke, area))
}

pub fn tri6_nodal_stress(
    xyz: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    plane_strain: bool,
) -> Result<Vec<[f64; 6]>> {
    let mut xy = [[0.0; 2]; 6];
    for i in 0..6 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let d = if plane_strain {
        d_plane_strain(e, nu)?
    } else {
        d_plane_stress(e, nu)?
    };
    let n = 12usize;
    let rst = [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.5, 0.0],
        [0.5, 0.5],
        [0.0, 0.5],
    ];
    let mut out = vec![[0.0; 6]; 6];
    for a in 0..6 {
        let (dndx, _, _) = tri6_dndx(&xy, rst[a][0], rst[a][1])?;
        let mut b = vec![0.0; 3 * n];
        fill_b2(&mut b, 6, &dndx);
        let s = sigma_from_b(&b, 3, n, &d, ue);
        out[a][0] = s[0];
        out[a][1] = s[1];
        out[a][3] = s[2];
    }
    Ok(out)
}

pub fn tri6_body_force(xy: &[[f64; 2]], bx: f64, by: f64, t: f64) -> Result<Vec<f64>> {
    let mut fe = vec![0.0; 12];
    let pts = [
        [1.0 / 6.0, 1.0 / 6.0],
        [2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0],
    ];
    let w = 1.0 / 6.0;
    for p in &pts {
        let (_, det, n) = tri6_dndx(xy, p[0], p[1])?;
        for a in 0..6 {
            fe[2 * a] += n[a] * bx * w * det * t;
            fe[2 * a + 1] += n[a] * by * w * det * t;
        }
    }
    Ok(fe)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quad8_partition_and_jac() {
        let xy = [
            [0.0, 0.0],
            [10.0, 0.0],
            [10.0, 10.0],
            [0.0, 10.0],
            [5.0, 0.0],
            [10.0, 5.0],
            [5.0, 10.0],
            [0.0, 5.0],
        ];
        let (n, _) = quad8_shape(0.0, 0.0);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "PU={sum}");
        let (_, det, _) = quad8_dndx(&xy, 0.0, 0.0).unwrap();
        assert!((det - 25.0).abs() < 1e-8, "det={det}");
        let (ke, area) = quad8_stiffness(&xy, 210000.0, 0.3, 1.0, false, false).unwrap();
        assert!((area - 100.0).abs() < 1e-6, "area={area}");
        let mut min_diag = f64::MAX;
        for i in 0..16 {
            min_diag = min_diag.min(ke[i * 16 + i]);
        }
        assert!(min_diag > 0.0, "min diag={min_diag}");
    }

    #[test]
    fn hex20_partition() {
        let (n, _) = hex20_shape(0.0, 0.0, 0.0);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "PU={sum}");
    }
}
