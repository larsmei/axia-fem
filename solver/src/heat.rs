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
        ElemKind::Quad4Ps | ElemKind::Quad4Pe => quad4_conductivity(xyz, k, area_or_th),
        _ => err(format!(
            "*HEAT TRANSFER: Element {} nicht implementiert (T3D2, B31, C3D8, CPS4).",
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
        ElemKind::Quad4Ps | ElemKind::Quad4Pe => quad4_body_heat(xyz, q, area_or_th),
        _ => Ok(vec![0.0; kind.nnodes()]),
    }
}
