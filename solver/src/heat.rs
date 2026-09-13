//! Steady (and lumped-capacity) heat conduction.
//!
//! Scalar DOF per node (NT). Conductivity from `*CONDUCTIVITY`.
//! Body/surface flux `*DFLUX`, convection `*FILM`, nodal `*CFLUX`.
//! Prescribed temperatures are `*BOUNDARY` DOF 11 (CalculiX NT).

use crate::elem::{hex8_dndx, invert2, G2, QUAD_XI};
use crate::error::{err, Result};
use crate::model::ElemKind;

const HEX_FACES: [[usize; 4]; 6] = [
    [0, 1, 2, 3],
    [4, 7, 6, 5],
    [0, 4, 5, 1],
    [1, 5, 6, 2],
    [2, 6, 7, 3],
    [3, 7, 4, 0],
];

fn hex_as_8(xyz: &[[f64; 3]]) -> Result<[[f64; 3]; 8]> {
    if xyz.len() < 8 {
        return err("C3D8 Wärmeleitung braucht 8 Knoten.");
    }
    let mut p = [[0.0; 3]; 8];
    for i in 0..8 {
        p[i] = xyz[i];
    }
    Ok(p)
}

/// Line element (T3D2 / B31): K = (k A / L) [1 -1; -1 1], lumped C = ρ cp A L / 2.
pub fn line_conductivity(xyz: &[[f64; 3]], k: f64, area: f64) -> Result<(Vec<f64>, f64)> {
    let nn = xyz.len().max(2);
    let i1 = if nn == 2 { 1 } else { nn - 1 };
    let dx = xyz[i1][0] - xyz[0][0];
    let dy = xyz[i1][1] - xyz[0][1];
    let dz = xyz[i1][2] - xyz[0][2];
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    if len < 1e-18 {
        return err("Wärmeleitung: Stablänge null.");
    }
    let coeff = k.max(0.0) * area.max(1e-30) / len;
    let mut ke = vec![0.0; nn * nn];
    ke[0] = coeff;
    ke[i1 * nn + i1] = coeff;
    ke[i1] = -coeff;
    ke[i1 * nn] = -coeff;
    Ok((ke, len * area))
}

pub fn line_body_heat(xyz: &[[f64; 3]], area: f64, q: f64) -> Vec<f64> {
    let nn = xyz.len().max(2);
    let i1 = if nn == 2 { 1 } else { nn - 1 };
    let dx = xyz[i1][0] - xyz[0][0];
    let dy = xyz[i1][1] - xyz[0][1];
    let dz = xyz[i1][2] - xyz[0][2];
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    let share = q * area * len / 2.0;
    let mut fe = vec![0.0; nn];
    fe[0] = share;
    fe[i1] = share;
    fe
}

pub fn hex8_conductivity(xyz: &[[f64; 3]], k: f64) -> Result<(Vec<f64>, f64)> {
    let p = hex_as_8(xyz)?;
    let mut ke = vec![0.0; 64];
    let mut vol = 0.0;
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            for &zeta in &pts {
                let (dndx, det, _) = hex8_dndx(&p, xi, eta, zeta)?;
                if det <= 0.0 {
                    return err("C3D8 Wärmeleitung: negative Jakobideterminante.");
                }
                vol += det;
                for a in 0..8 {
                    for b in 0..8 {
                        let dot = dndx[a][0] * dndx[b][0]
                            + dndx[a][1] * dndx[b][1]
                            + dndx[a][2] * dndx[b][2];
                        ke[a * 8 + b] += k * det * dot;
                    }
                }
            }
        }
    }
    Ok((ke, vol))
}

pub fn hex8_body_heat(xyz: &[[f64; 3]], q: f64) -> Result<Vec<f64>> {
    let p = hex_as_8(xyz)?;
    let mut fe = vec![0.0; 8];
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            for &zeta in &pts {
                let (_, det, n) = hex8_dndx(&p, xi, eta, zeta)?;
                if det <= 0.0 {
                    continue;
                }
                for a in 0..8 {
                    fe[a] += q * n[a] * det;
                }
            }
        }
    }
    Ok(fe)
}

fn face_shape(xi: f64, eta: f64) -> ([f64; 4], [f64; 4], [f64; 4]) {
    let mut n = [0.0; 4];
    let mut dnxi = [0.0; 4];
    let mut dneta = [0.0; 4];
    for i in 0..4 {
        let x = QUAD_XI[i][0];
        let e = QUAD_XI[i][1];
        n[i] = 0.25 * (1.0 + x * xi) * (1.0 + e * eta);
        dnxi[i] = 0.25 * x * (1.0 + e * eta);
        dneta[i] = 0.25 * e * (1.0 + x * xi);
    }
    (n, dnxi, dneta)
}

/// Face flux + optional film. `ke` is 8×8 added convection, `fe` is 8 nodal heat.
pub fn hex8_face_heat(
    xyz: &[[f64; 3]],
    face: i32,
    flux: f64,
    h: f64,
    t_inf: f64,
) -> Result<(Vec<f64>, Vec<f64>)> {
    if !(1..=6).contains(&face) {
        return err(format!("Ungültige C3D8-Fläche S{face}"));
    }
    let p = hex_as_8(xyz)?;
    let fi = (face - 1) as usize;
    let mut fe = vec![0.0; 8];
    let mut ke = vec![0.0; 64];
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            let (n, dnxi, dneta) = face_shape(xi, eta);
            let mut rxi = [0.0; 3];
            let mut reta = [0.0; 3];
            for a in 0..4 {
                let q = p[HEX_FACES[fi][a]];
                for k in 0..3 {
                    rxi[k] += dnxi[a] * q[k];
                    reta[k] += dneta[a] * q[k];
                }
            }
            let nx = rxi[1] * reta[2] - rxi[2] * reta[1];
            let ny = rxi[2] * reta[0] - rxi[0] * reta[2];
            let nz = rxi[0] * reta[1] - rxi[1] * reta[0];
            let jac = (nx * nx + ny * ny + nz * nz).sqrt();
            let qn = flux + h * t_inf;
            for a in 0..4 {
                let ia = HEX_FACES[fi][a];
                fe[ia] += n[a] * qn * jac;
                if h.abs() > 0.0 {
                    for b in 0..4 {
                        let ib = HEX_FACES[fi][b];
                        ke[ia * 8 + ib] += h * n[a] * n[b] * jac;
                    }
                }
            }
        }
    }
    Ok((ke, fe))
}

pub fn quad4_conductivity(xyz: &[[f64; 3]], k: f64, th: f64) -> Result<(Vec<f64>, f64)> {
    let mut xy = [[0.0; 2]; 4];
    for i in 0..4 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let mut ke = vec![0.0; 16];
    let mut area = 0.0;
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            let mut n = [0.0; 4];
            let mut dn = [[0.0; 2]; 4];
            for i in 0..4 {
                let x = QUAD_XI[i][0];
                let e = QUAD_XI[i][1];
                n[i] = 0.25 * (1.0 + x * xi) * (1.0 + e * eta);
                dn[i][0] = 0.25 * x * (1.0 + e * eta);
                dn[i][1] = 0.25 * e * (1.0 + x * xi);
            }
            let mut j = [[0.0; 2]; 2];
            for a in 0..4 {
                j[0][0] += dn[a][0] * xy[a][0];
                j[0][1] += dn[a][1] * xy[a][0];
                j[1][0] += dn[a][0] * xy[a][1];
                j[1][1] += dn[a][1] * xy[a][1];
            }
            let (inv, det) = invert2(j)?;
            if det <= 0.0 {
                return err("CPS4 Wärmeleitung: negative Jakobideterminante.");
            }
            area += det;
            let mut dndx = [[0.0; 2]; 4];
            for a in 0..4 {
                dndx[a][0] = inv[0][0] * dn[a][0] + inv[1][0] * dn[a][1];
                dndx[a][1] = inv[0][1] * dn[a][0] + inv[1][1] * dn[a][1];
            }
            let w = k * th * det;
            for a in 0..4 {
                for b in 0..4 {
                    ke[a * 4 + b] += w * (dndx[a][0] * dndx[b][0] + dndx[a][1] * dndx[b][1]);
                }
            }
            let _ = n;
        }
    }
    Ok((ke, area * th))
}

pub fn quad4_body_heat(xyz: &[[f64; 3]], q: f64, th: f64) -> Result<Vec<f64>> {
    let mut xy = [[0.0; 2]; 4];
    for i in 0..4 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let mut fe = vec![0.0; 4];
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            let mut n = [0.0; 4];
            let mut dn = [[0.0; 2]; 4];
            for i in 0..4 {
                let x = QUAD_XI[i][0];
                let e = QUAD_XI[i][1];
                n[i] = 0.25 * (1.0 + x * xi) * (1.0 + e * eta);
                dn[i][0] = 0.25 * x * (1.0 + e * eta);
                dn[i][1] = 0.25 * e * (1.0 + x * xi);
            }
            let mut j = [[0.0; 2]; 2];
            for a in 0..4 {
                j[0][0] += dn[a][0] * xy[a][0];
                j[0][1] += dn[a][1] * xy[a][0];
                j[1][0] += dn[a][0] * xy[a][1];
                j[1][1] += dn[a][1] * xy[a][1];
            }
            let (_, det) = invert2(j)?;
            if det <= 0.0 {
                continue;
            }
            for a in 0..4 {
                fe[a] += q * th * n[a] * det;
            }
        }
    }
    Ok(fe)
}

fn hex20_conductivity(xyz: &[[f64; 3]], k: f64, reduced: bool) -> Result<(Vec<f64>, f64)> {
    let nn = 20usize;
    let mut ke = vec![0.0; nn * nn];
    let mut vol = 0.0;
    for (xi, eta, zeta, w0) in crate::quadratic::hex_gauss(reduced) {
        let (dndx, det, _) = crate::quadratic::hex20_dndx(xyz, xi, eta, zeta)?;
        if det <= 0.0 {
            return err("C3D20 Wärmeleitung: negative Jakobideterminante.");
        }
        let w = k * w0 * det;
        vol += w0 * det;
        for a in 0..nn {
            for b in 0..nn {
                ke[a * nn + b] += w
                    * (dndx[a][0] * dndx[b][0] + dndx[a][1] * dndx[b][1] + dndx[a][2] * dndx[b][2]);
            }
        }
    }
    Ok((ke, vol))
}

fn hex20_body_heat(xyz: &[[f64; 3]], q: f64, reduced: bool) -> Result<Vec<f64>> {
    let mut fe = vec![0.0; 20];
    for (xi, eta, zeta, w0) in crate::quadratic::hex_gauss(reduced) {
        let (_, det, n) = crate::quadratic::hex20_dndx(xyz, xi, eta, zeta)?;
        if det <= 0.0 {
            continue;
        }
        for a in 0..20 {
            fe[a] += q * n[a] * w0 * det;
        }
    }
    Ok(fe)
}

fn tet4_conductivity(xyz: &[[f64; 3]], k: f64) -> Result<(Vec<f64>, f64)> {
    let mut j = [[0.0; 3]; 3];
    for p in 0..3 {
        for q in 0..3 {
            j[q][p] = xyz[p + 1][q] - xyz[0][q];
        }
    }
    let (inv, det) = crate::elem::invert3(j)?;
    if det <= 0.0 {
        return err("C3D4 Wärmeleitung: negative Jakobideterminante.");
    }
    let vol = det / 6.0;
    let dn = [
        [-1.0, -1.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mut dndx = [[0.0; 3]; 4];
    for a in 0..4 {
        for i in 0..3 {
            dndx[a][i] = inv[0][i] * dn[a][0] + inv[1][i] * dn[a][1] + inv[2][i] * dn[a][2];
        }
    }
    let mut ke = vec![0.0; 16];
    for a in 0..4 {
        for b in 0..4 {
            ke[a * 4 + b] = k
                * vol
                * (dndx[a][0] * dndx[b][0] + dndx[a][1] * dndx[b][1] + dndx[a][2] * dndx[b][2]);
        }
    }
    Ok((ke, vol))
}

fn tet10_conductivity(xyz: &[[f64; 3]], k: f64) -> Result<(Vec<f64>, f64)> {
    let nn = 10usize;
    let mut ke = vec![0.0; nn * nn];
    let mut vol = 0.0;
    let a = 0.5854101966249685;
    let b = 0.1381966011250105;
    let w0 = 1.0 / 24.0;
    let pts = [[b, b, b], [a, b, b], [b, a, b], [b, b, a]];
    for p in &pts {
        let (dndx, det, _) = crate::quadratic::tet10_dndx(xyz, p[0], p[1], p[2])?;
        if det <= 0.0 {
            return err("C3D10 Wärmeleitung: negative Jakobideterminante.");
        }
        let w = k * w0 * det;
        vol += w0 * det;
        for i in 0..nn {
            for j in 0..nn {
                ke[i * nn + j] += w
                    * (dndx[i][0] * dndx[j][0] + dndx[i][1] * dndx[j][1] + dndx[i][2] * dndx[j][2]);
            }
        }
    }
    Ok((ke, vol))
}

fn wedge6_conductivity(xyz: &[[f64; 3]], k: f64) -> Result<(Vec<f64>, f64)> {
    fn shape(xi: f64, eta: f64, zeta: f64) -> ([f64; 6], [[f64; 3]; 6]) {
        let l1 = xi;
        let l2 = eta;
        let l3 = 1.0 - xi - eta;
        let zm = 0.5 * (1.0 - zeta);
        let zp = 0.5 * (1.0 + zeta);
        let n = [l1 * zm, l2 * zm, l3 * zm, l1 * zp, l2 * zp, l3 * zp];
        let mut dn = [[0.0; 3]; 6];
        dn[0] = [zm, 0.0, -0.5 * l1];
        dn[1] = [0.0, zm, -0.5 * l2];
        dn[2] = [-zm, -zm, -0.5 * l3];
        dn[3] = [zp, 0.0, 0.5 * l1];
        dn[4] = [0.0, zp, 0.5 * l2];
        dn[5] = [-zp, -zp, 0.5 * l3];
        (n, dn)
    }
    let nn = 6usize;
    let mut ke = vec![0.0; nn * nn];
    let mut vol = 0.0;
    let tri = [
        [1.0 / 6.0, 1.0 / 6.0],
        [2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0],
    ];
    let wtri = 1.0 / 6.0;
    for t in &tri {
        for &zeta in &[-G2, G2] {
            let (n, dn) = shape(t[0], t[1], zeta);
            let _ = n;
            let mut j = [[0.0; 3]; 3];
            for a in 0..6 {
                for p in 0..3 {
                    for q in 0..3 {
                        j[q][p] += dn[a][p] * xyz[a][q];
                    }
                }
            }
            let (inv, det) = crate::elem::invert3(j)?;
            if det <= 0.0 {
                return err("C3D6 Wärmeleitung: negative Jakobideterminante.");
            }
            let w = k * det * wtri;
            vol += det * wtri;
            let mut dndx = [[0.0; 3]; 6];
            for a in 0..6 {
                for i in 0..3 {
                    dndx[a][i] = inv[0][i] * dn[a][0] + inv[1][i] * dn[a][1] + inv[2][i] * dn[a][2];
                }
            }
            for a in 0..6 {
                for b in 0..6 {
                    ke[a * nn + b] += w
                        * (dndx[a][0] * dndx[b][0]
                            + dndx[a][1] * dndx[b][1]
                            + dndx[a][2] * dndx[b][2]);
                }
            }
        }
    }
    Ok((ke, vol))
}

fn wedge15_conductivity(xyz: &[[f64; 3]], k: f64) -> Result<(Vec<f64>, f64)> {
    let nn = 15usize;
    let mut ke = vec![0.0; nn * nn];
    let mut vol = 0.0;
    let tri = [
        [1.0 / 6.0, 1.0 / 6.0],
        [2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0],
    ];
    let wtri = 1.0 / 6.0;
    const G3: [f64; 3] = [-0.7745966692414834, 0.0, 0.7745966692414834];
    const W3: [f64; 3] = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
    for t in &tri {
        for i in 0..3 {
            let (dndx, det, _) = crate::extra::wedge15_dndx(xyz, t[0], t[1], G3[i])?;
            if det <= 0.0 {
                return err("C3D15 Wärmeleitung: negative Jakobideterminante.");
            }
            let w = k * det * wtri * W3[i];
            vol += det * wtri * W3[i];
            for a in 0..nn {
                for b in 0..nn {
                    ke[a * nn + b] += w
                        * (dndx[a][0] * dndx[b][0]
                            + dndx[a][1] * dndx[b][1]
                            + dndx[a][2] * dndx[b][2]);
                }
            }
        }
    }
    Ok((ke, vol))
}

fn quad8_conductivity(xyz: &[[f64; 3]], k: f64, th: f64, reduced: bool) -> Result<(Vec<f64>, f64)> {
    let mut xy = [[0.0; 2]; 8];
    for i in 0..8 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let nn = 8usize;
    let mut ke = vec![0.0; nn * nn];
    let mut area = 0.0;
    let gps: Vec<(f64, f64, f64)> = if reduced {
        let mut o = Vec::new();
        for &xi in &[-G2, G2] {
            for &eta in &[-G2, G2] {
                o.push((xi, eta, 1.0));
            }
        }
        o
    } else {
        const G3: [f64; 3] = [-0.7745966692414834, 0.0, 0.7745966692414834];
        const W3: [f64; 3] = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
        let mut o = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                o.push((G3[i], G3[j], W3[i] * W3[j]));
            }
        }
        o
    };
    for (xi, eta, w0) in gps {
        let (dndx, det, _) = crate::quadratic::quad8_dndx(&xy, xi, eta)?;
        if det <= 0.0 {
            return err("CPS8 Wärmeleitung: negative Jakobideterminante.");
        }
        let w = k * th * w0 * det;
        area += w0 * det;
        for a in 0..nn {
            for b in 0..nn {
                ke[a * nn + b] += w * (dndx[a][0] * dndx[b][0] + dndx[a][1] * dndx[b][1]);
            }
        }
    }
    Ok((ke, area * th))
}

fn shell_quad_conductivity(xyz: &[[f64; 3]], k: f64, th: f64) -> Result<(Vec<f64>, f64)> {
    // in-plane conduction in the shell tangent plane
    let (e1, e2, _) = crate::shell::local_frame(xyz, 4)?;
    let mut xy = [[0.0; 2]; 4];
    for i in 0..4 {
        let d = [
            xyz[i][0] - xyz[0][0],
            xyz[i][1] - xyz[0][1],
            xyz[i][2] - xyz[0][2],
        ];
        xy[i] = [
            d[0] * e1[0] + d[1] * e1[1] + d[2] * e1[2],
            d[0] * e2[0] + d[1] * e2[1] + d[2] * e2[2],
        ];
    }
    let mut xyz2 = vec![[0.0; 3]; 4];
    for i in 0..4 {
        xyz2[i] = [xy[i][0], xy[i][1], 0.0];
    }
    quad4_conductivity(&xyz2, k, th)
}

fn tri3_conductivity(xyz: &[[f64; 3]], k: f64, th: f64) -> Result<(Vec<f64>, f64)> {
    if xyz.len() < 3 {
        return err("Dreieck-Wärmeleitung braucht 3 Knoten.");
    }
    let two_a = xyz[0][0] * (xyz[1][1] - xyz[2][1])
        + xyz[1][0] * (xyz[2][1] - xyz[0][1])
        + xyz[2][0] * (xyz[0][1] - xyz[1][1]);
    if two_a.abs() < 1e-18 {
        return err("Dreieck-Wärmeleitung: Fläche null.");
    }
    let area = two_a.abs() / 2.0;
    let dndx = [
        [(xyz[1][1] - xyz[2][1]) / two_a, (xyz[2][0] - xyz[1][0]) / two_a],
        [(xyz[2][1] - xyz[0][1]) / two_a, (xyz[0][0] - xyz[2][0]) / two_a],
        [(xyz[0][1] - xyz[1][1]) / two_a, (xyz[1][0] - xyz[0][0]) / two_a],
    ];
    let mut ke = vec![0.0; 9];
    let w = k.max(0.0) * th.max(1e-30) * area;
    for a in 0..3 {
        for b in 0..3 {
            ke[a * 3 + b] = w * (dndx[a][0] * dndx[b][0] + dndx[a][1] * dndx[b][1]);
        }
    }
    Ok((ke, area * th))
}

fn tri6_conductivity(xyz: &[[f64; 3]], k: f64, th: f64) -> Result<(Vec<f64>, f64)> {
    if xyz.len() < 6 {
        return err("CPS6-Wärmeleitung braucht 6 Knoten.");
    }
    let mut xy = [[0.0; 2]; 6];
    for i in 0..6 {
        xy[i] = [xyz[i][0], xyz[i][1]];
    }
    let nn = 6usize;
    let mut ke = vec![0.0; nn * nn];
    let mut area = 0.0;
    let pts = [
        [1.0 / 6.0, 1.0 / 6.0],
        [2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0],
    ];
    let w0 = 1.0 / 6.0;
    for p in &pts {
        let (dndx, det, _) = crate::quadratic::tri6_dndx(&xy, p[0], p[1])?;
        if det <= 0.0 {
            return err("CPS6 Wärmeleitung: negative Jakobideterminante.");
        }
        let w = k.max(0.0) * th.max(1e-30) * w0 * det;
        area += w0 * det;
        for a in 0..nn {
            for b in 0..nn {
                ke[a * nn + b] += w * (dndx[a][0] * dndx[b][0] + dndx[a][1] * dndx[b][1]);
            }
        }
    }
    Ok((ke, area * th))
}

pub fn element_conductivity(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    k: f64,
    area_or_th: f64,
) -> Result<(Vec<f64>, f64)> {
    match kind {
        ElemKind::Truss2 | ElemKind::Truss3 | ElemKind::Beam31 | ElemKind::Beam32 => {
            line_conductivity(xyz, k, area_or_th)
        }
        ElemKind::Hex8 | ElemKind::Hex8I | ElemKind::Hex8R => hex8_conductivity(xyz, k),
        ElemKind::Hex20 | ElemKind::Hex20R => hex20_conductivity(xyz, k, kind.reduced_int()),
        ElemKind::Tet4 => tet4_conductivity(xyz, k),
        ElemKind::Tet10 | ElemKind::Tet10T => tet10_conductivity(xyz, k),
        ElemKind::Wedge6 => wedge6_conductivity(xyz, k),
        ElemKind::Wedge15 => wedge15_conductivity(xyz, k),
        ElemKind::Quad4Ps | ElemKind::Quad4Pe => quad4_conductivity(xyz, k, area_or_th),
        ElemKind::Quad8Ps | ElemKind::Quad8Pe | ElemKind::Quad8RPs | ElemKind::Quad8RPe => {
            quad8_conductivity(xyz, k, area_or_th, kind.reduced_int())
        }
        ElemKind::Shell4 | ElemKind::Shell4R | ElemKind::Mem4 | ElemKind::Mem4R => {
            shell_quad_conductivity(xyz, k, area_or_th)
        }
        ElemKind::Cax4 | ElemKind::Cax4R => {
            let r = xyz.iter().map(|p| p[0]).sum::<f64>() / xyz.len().max(1) as f64;
            shell_quad_conductivity(xyz, k, 2.0 * std::f64::consts::PI * r.max(1e-16))
        }
        ElemKind::Cax8 | ElemKind::Cax8R | ElemKind::Shell8 | ElemKind::Shell8R | ElemKind::Mem8 => {
            let r = xyz.iter().map(|p| p[0]).sum::<f64>() / xyz.len().max(1) as f64;
            let th = if kind.is_axisym() {
                2.0 * std::f64::consts::PI * r.max(1e-16)
            } else {
                area_or_th
            };
            quad8_conductivity(xyz, k, th, kind.reduced_int())
        }
        ElemKind::Cax3 | ElemKind::Tri3Ps | ElemKind::Tri3Pe | ElemKind::Shell3 | ElemKind::Mem3 => {
            let th = if kind.is_axisym() {
                let r = xyz.iter().map(|p| p[0]).sum::<f64>() / xyz.len().max(1) as f64;
                2.0 * std::f64::consts::PI * r.max(1e-16)
            } else {
                area_or_th
            };
            tri3_conductivity(xyz, k, th)
        }
        ElemKind::Cax6 | ElemKind::Tri6Ps | ElemKind::Tri6Pe | ElemKind::Shell6 | ElemKind::Mem6 => {
            let th = if kind.is_axisym() {
                let r = xyz.iter().map(|p| p[0]).sum::<f64>() / xyz.len().max(1) as f64;
                2.0 * std::f64::consts::PI * r.max(1e-16)
            } else {
                area_or_th
            };
            tri6_conductivity(xyz, k, th)
        }
        ElemKind::Mass | ElemKind::RotaryI | ElemKind::DashpotA | ElemKind::GapUni => {
            Ok((vec![0.0; kind.nnodes() * kind.nnodes()], 0.0))
        }
        _ => err(format!(
            "*HEAT TRANSFER: Element {} nicht implementiert.",
            kind.ccx_name()
        )),
    }
}

pub fn element_body_heat(kind: ElemKind, xyz: &[[f64; 3]], q: f64, area_or_th: f64) -> Result<Vec<f64>> {
    match kind {
        ElemKind::Truss2 | ElemKind::Truss3 | ElemKind::Beam31 | ElemKind::Beam32 => {
            Ok(line_body_heat(xyz, area_or_th, q))
        }
        ElemKind::Hex8 | ElemKind::Hex8I | ElemKind::Hex8R => hex8_body_heat(xyz, q),
        ElemKind::Hex20 | ElemKind::Hex20R => hex20_body_heat(xyz, q, kind.reduced_int()),
        ElemKind::Quad4Ps | ElemKind::Quad4Pe => quad4_body_heat(xyz, q, area_or_th),
        _ => Ok(vec![0.0; kind.nnodes()]),
    }
}
