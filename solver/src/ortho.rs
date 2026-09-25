//! Orthotropic and anisotropic plate/solid laws for MYSTRAN MAT2, MAT8, MAT9 and PCOMP.
//!
//! Plane-stress Voigt order is [εx, εy, γxy]. The 3D MAT9 order matches the manual
//! and Axia: [εx, εy, εz, γxy, γyz, γzx].

use crate::error::{err, Result};

const K_SHEAR: f64 = 5.0 / 6.0;

/// Membrane A, bending D and transverse shear, already in the element frame.
/// Coupling B is not stored: MITC4 has no ABD coupling.
#[derive(Clone, Debug)]
pub struct ShellLaw {
    pub a: [f64; 9],
    pub d: [f64; 9],
    pub ds: [f64; 4],
    /// `E·h` for the drilling stabilizer. MITC4 multiplies its own drill factor.
    pub drill_eh: f64,
}

/// Isotropic plane-stress plate. `bi` is PSHELL 12I/T³.
pub fn isotropic_shell(e: f64, nu: f64, h: f64, bi: f64) -> Result<ShellLaw> {
    if nu >= 1.0 || nu <= -1.0 {
        return err(format!("Ungültige Querkontraktion nu={nu}"));
    }
    if e <= 0.0 {
        return err(format!("Ungültiger E-Modul E={e}"));
    }
    if h <= 0.0 {
        return err("SHELL SECTION: Dicke muss positiv sein.");
    }
    let c = e / (1.0 - nu * nu);
    let mut q = [0.0; 9];
    q[0] = c;
    q[1] = c * nu;
    q[3] = c * nu;
    q[4] = c;
    q[8] = c * (1.0 - nu) / 2.0;
    let gsh = e / (2.0 * (1.0 + nu)) * K_SHEAR * h;
    Ok(homogenous_q(&q, h, bi, [gsh, 0.0, 0.0, gsh], e))
}

/// Orthotropic plane-stress Q in the material axes. `nu12` is −ε2/ε1 under σ1.
pub fn q_ortho(e1: f64, e2: f64, nu12: f64, g12: f64) -> Result<[f64; 9]> {
    if e1 <= 0.0 || e2 <= 0.0 {
        return err(format!("MAT8: E1={e1} und E2={e2} müssen positiv sein."));
    }
    if g12 < 0.0 {
        return err(format!("MAT8: G12={g12} darf nicht negativ sein."));
    }
    let nu21 = nu12 * e2 / e1;
    let den = 1.0 - nu12 * nu21;
    if den <= 1e-12 {
        return err(format!("MAT8: 1−ν12·ν21={den} ist nicht positiv."));
    }
    let mut q = [0.0; 9];
    q[0] = e1 / den;
    q[1] = nu12 * e2 / den;
    q[3] = q[1];
    q[4] = e2 / den;
    q[8] = g12;
    Ok(q)
}

/// Rotate a plane-stress matrix by `theta_deg` (material 1-axis from element x).
pub fn rotate_q(q: &[f64; 9], theta_deg: f64) -> [f64; 9] {
    let t = theta_deg.to_radians();
    let c = t.cos();
    let s = t.sin();
    let c2 = c * c;
    let s2 = s * s;
    let cs = c * s;
    let q11 = q[0];
    let q12 = q[1];
    let q16 = q[2];
    let q22 = q[4];
    let q26 = q[5];
    let q66 = q[8];
    // Standard Q̄. At 90° this swaps Q11 and Q22 when Q16 = Q26 = 0.
    let q11b = q11 * c2 * c2
        + 2.0 * (q12 + 2.0 * q66) * s2 * c2
        + q22 * s2 * s2
        + 4.0 * q16 * c2 * cs
        + 4.0 * q26 * s2 * cs;
    let q22b = q11 * s2 * s2
        + 2.0 * (q12 + 2.0 * q66) * s2 * c2
        + q22 * c2 * c2
        - 4.0 * q16 * s2 * cs
        - 4.0 * q26 * c2 * cs;
    let q12b = (q11 + q22 - 4.0 * q66) * s2 * c2
        + q12 * (s2 * s2 + c2 * c2)
        + 2.0 * (q26 - q16) * cs * (c2 - s2);
    let q16b = (q11 - q12 - 2.0 * q66) * cs * c2
        + (q12 - q22 + 2.0 * q66) * cs * s2
        + q16 * c2 * (c2 - 3.0 * s2)
        + q26 * s2 * (3.0 * c2 - s2);
    let q26b = (q11 - q12 - 2.0 * q66) * cs * s2
        + (q12 - q22 + 2.0 * q66) * cs * c2
        + q16 * s2 * (s2 - 3.0 * c2)
        + q26 * c2 * (c2 - 3.0 * s2);
    let q66b = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * s2 * c2
        + q66 * (s2 * s2 + c2 * c2)
        + 2.0 * (q16 - q26) * cs * (c2 - s2);
    let mut o = [0.0; 9];
    o[0] = q11b;
    o[1] = q12b;
    o[2] = q16b;
    o[3] = q12b;
    o[4] = q22b;
    o[5] = q26b;
    o[6] = q16b;
    o[7] = q26b;
    o[8] = q66b;
    o
}

/// Rotate transverse-shear moduli (G1z along the material 1-axis).
pub fn rotate_shear(g1: f64, g2: f64, theta_deg: f64) -> [f64; 4] {
    let t = theta_deg.to_radians();
    let c = t.cos();
    let s = t.sin();
    let c2 = c * c;
    let s2 = s * s;
    let cs = c * s;
    let xx = g1 * c2 + g2 * s2;
    let yy = g1 * s2 + g2 * c2;
    let xy = (g1 - g2) * cs;
    [xx, xy, xy, yy]
}

fn homogenous_q(q: &[f64; 9], h: f64, bi: f64, ds: [f64; 4], drill_e: f64) -> ShellLaw {
    let mut a = [0.0; 9];
    let mut d = [0.0; 9];
    for i in 0..9 {
        a[i] = q[i] * h;
        d[i] = q[i] * h * h * h / 12.0 * bi;
    }
    ShellLaw {
        a,
        d,
        ds,
        drill_eh: drill_e * h,
    }
}

/// Homogeneous MAT2/MAT8 plate. `g1z`/`g2z` ≤ 0 means zero shear flexibility.
pub fn homogeneous_plate(
    q: &[f64; 9],
    h: f64,
    bi: f64,
    theta_deg: f64,
    g1z: f64,
    g2z: f64,
    drill_e: f64,
) -> Result<ShellLaw> {
    if h <= 0.0 {
        return err("Plattendicke muss positiv sein.");
    }
    let qb = rotate_q(q, theta_deg);
    let ds = transverse_shear(g1z, g2z, theta_deg, h, drill_e);
    Ok(homogenous_q(&qb, h, bi, ds, drill_e))
}

pub fn transverse_shear(g1z: f64, g2z: f64, theta_deg: f64, h: f64, drill_e: f64) -> [f64; 4] {
    let big = drill_e.abs().max(1.0) * h.abs().max(1e-12) * 1.0e6;
    let g1 = if g1z > 0.0 { g1z * K_SHEAR * h } else { big };
    let g2 = if g2z > 0.0 { g2z * K_SHEAR * h } else { big };
    rotate_shear(g1, g2, theta_deg)
}

#[derive(Clone, Copy, Debug)]
pub struct PlyQ {
    pub q: [f64; 9],
    pub t: f64,
    pub theta_deg: f64,
    pub g1z: f64,
    pub g2z: f64,
    pub drill_e: f64,
}

/// ABD from plies already expanded (SYM applied) and measured from `z0`.
/// Returns the law and whether |B| is large enough to mention.
pub fn laminate(plies: &[PlyQ], z0: f64) -> Result<(ShellLaw, bool)> {
    if plies.is_empty() {
        return err("PCOMP ohne Lagen.");
    }
    let mut a = [0.0; 9];
    let mut b = [0.0; 9];
    let mut d = [0.0; 9];
    let mut ds = [0.0; 4];
    let mut z = z0;
    let mut drill_eh = 0.0;
    let mut thick = 0.0;
    for p in plies {
        if p.t <= 0.0 {
            return err("PCOMP: Lagendicke muss positiv sein.");
        }
        let z1 = z;
        let z2 = z + p.t;
        let q = rotate_q(&p.q, p.theta_deg);
        let dz = p.t;
        let dz2 = 0.5 * (z2 * z2 - z1 * z1);
        let dz3 = (z2 * z2 * z2 - z1 * z1 * z1) / 3.0;
        for i in 0..9 {
            a[i] += q[i] * dz;
            b[i] += q[i] * dz2;
            d[i] += q[i] * dz3;
        }
        let s = transverse_shear(p.g1z, p.g2z, p.theta_deg, p.t, p.drill_e);
        for i in 0..4 {
            ds[i] += s[i];
        }
        drill_eh += p.drill_e * p.t;
        thick += p.t;
        z = z2;
    }
    if thick <= 0.0 {
        return err("PCOMP: Gesamtdicke muss positiv sein.");
    }
    let bn = frob9(&b);
    let scale = frob9(&a).sqrt() * frob9(&d).sqrt();
    let coupled = bn > 1e-8 * scale.max(1.0);
    Ok((
        ShellLaw {
            a,
            d,
            ds,
            drill_eh,
        },
        coupled,
    ))
}

fn frob9(m: &[f64; 9]) -> f64 {
    m.iter().map(|v| v * v).sum::<f64>().sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn theta_90_swaps_e1_and_e2() {
        let q = q_ortho(2.0, 3.0, 0.0, 1.0).unwrap();
        let qb = rotate_q(&q, 90.0);
        assert!((qb[0] - 3.0).abs() < 1e-12, "Q11 {}", qb[0]);
        assert!((qb[4] - 2.0).abs() < 1e-12, "Q22 {}", qb[4]);
        assert!(qb[1].abs() < 1e-12);
        assert!(qb[2].abs() < 1e-12);
    }
}
