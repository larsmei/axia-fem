//! Small-strain J2 plasticity with isotropic hardening for continuum (C3D*).
//! CalculiX `*PLASTIC` curve (σy, ε̄p) is stored as (peeq, sy).

use crate::elem::{d_iso_3d, fill_b3, gemm_bt_d_b, hex8_dndx, invert3, G2};
use crate::error::{err, Result};
use crate::material::yield_from_curve;
use crate::model::ElemKind;
use crate::quadratic;

#[derive(Clone, Copy, Debug)]
pub struct GpHist {
    pub ep: [f64; 6],
    pub alpha: f64,
}

impl Default for GpHist {
    fn default() -> Self {
        Self {
            ep: [0.0; 6],
            alpha: 0.0,
        }
    }
}

pub struct PlasticElem {
    pub ke: Vec<f64>,
    pub fe: Vec<f64>,
    pub stress: [f64; 6],
    pub strain: [f64; 6],
    pub peeq: f64,
    #[allow(dead_code)]
    pub ngp: usize,
}

fn matvec6(c: &[f64], e: &[f64; 6]) -> [f64; 6] {
    let mut s = [0.0; 6];
    for i in 0..6 {
        let mut v = 0.0;
        for j in 0..6 {
            v += c[i * 6 + j] * e[j];
        }
        s[i] = v;
    }
    s
}

fn vm_q(s: &[f64; 6]) -> f64 {
    let p = (s[0] + s[1] + s[2]) / 3.0;
    let d0 = s[0] - p;
    let d1 = s[1] - p;
    let d2 = s[2] - p;
    let j2 = 0.5 * (d0 * d0 + d1 * d1 + d2 * d2) + s[3] * s[3] + s[4] * s[4] + s[5] * s[5];
    (3.0 * j2.max(0.0)).sqrt()
}

fn deviator(s: &[f64; 6]) -> ([f64; 6], f64) {
    let p = (s[0] + s[1] + s[2]) / 3.0;
    (
        [s[0] - p, s[1] - p, s[2] - p, s[3], s[4], s[5]],
        p,
    )
}

fn snorm(dev: &[f64; 6]) -> f64 {
    (dev[0] * dev[0]
        + dev[1] * dev[1]
        + dev[2] * dev[2]
        + 2.0 * (dev[3] * dev[3] + dev[4] * dev[4] + dev[5] * dev[5]))
        .sqrt()
}

/// Radial return. `eps` and `ep` are engineering Voigt. Returns (σ, Cep, ep_new, α_new).
pub fn j2_return(
    eps: &[f64; 6],
    hist: &GpHist,
    e: f64,
    nu: f64,
    curve: &[(f64, f64)],
) -> Result<([f64; 6], [f64; 36], GpHist)> {
    let c = d_iso_3d(e, nu)?;
    let mut ee = [0.0; 6];
    for i in 0..6 {
        ee[i] = eps[i] - hist.ep[i];
    }
    let sig_tr = matvec6(&c, &ee);
    if curve.is_empty() {
        return Ok((sig_tr, c, *hist));
    }
    let q_tr = vm_q(&sig_tr);
    let (sy0, h0) = yield_from_curve(curve, hist.alpha);
    if q_tr <= sy0 * (1.0 + 1e-12) + 1e-10 {
        return Ok((sig_tr, c, *hist));
    }
    let g = e / (2.0 * (1.0 + nu));
    let kbulk = e / (3.0 * (1.0 - 2.0 * nu));
    let mut dg = ((q_tr - sy0) / (3.0 * g + h0.max(0.0))).max(0.0);
    let mut h = h0;
    let mut sy = sy0;
    for _ in 0..40 {
        let yh = yield_from_curve(curve, hist.alpha + dg);
        sy = yh.0;
        h = yh.1;
        let r = q_tr - 3.0 * g * dg - sy;
        if r.abs() < 1e-12 * (1.0 + q_tr) {
            break;
        }
        dg = (dg + r / (3.0 * g + h.max(0.0))).max(0.0);
    }
    let (dev_tr, p) = deviator(&sig_tr);
    let scale = if q_tr > 1e-16 {
        1.0 - 3.0 * g * dg / q_tr
    } else {
        0.0
    };
    let mut sig = [
        scale * dev_tr[0] + p,
        scale * dev_tr[1] + p,
        scale * dev_tr[2] + p,
        scale * dev_tr[3],
        scale * dev_tr[4],
        scale * dev_tr[5],
    ];
    let mut ep = hist.ep;
    if q_tr > 1e-16 {
        // Δεp_ij = Δγ (3/2) s_ij / q   (γp_xy = 2 Δεp_xy)
        let f = dg * 1.5 / q_tr;
        ep[0] += f * dev_tr[0];
        ep[1] += f * dev_tr[1];
        ep[2] += f * dev_tr[2];
        ep[3] += 2.0 * f * dev_tr[3];
        ep[4] += 2.0 * f * dev_tr[4];
        ep[5] += 2.0 * f * dev_tr[5];
    }
    let nrm = snorm(&dev_tr).max(1e-30);
    let n = [
        dev_tr[0] / nrm,
        dev_tr[1] / nrm,
        dev_tr[2] / nrm,
        dev_tr[3] / nrm,
        dev_tr[4] / nrm,
        dev_tr[5] / nrm,
    ];
    let theta = 1.0 - 3.0 * g * dg / q_tr.max(1e-30);
    let theta_bar = 3.0 * g / (3.0 * g + h.max(0.0)) - (1.0 - theta);
    let mut cep = [0.0; 36];
    for i in 0..3 {
        for j in 0..3 {
            cep[i * 6 + j] = kbulk;
        }
    }
    let twog = 2.0 * g * theta;
    for i in 0..3 {
        for j in 0..3 {
            cep[i * 6 + j] += twog * (if i == j { 1.0 } else { 0.0 } - 1.0 / 3.0);
        }
    }
    for i in 3..6 {
        cep[i * 6 + i] += twog * 0.5;
    }
    let fac = 2.0 * g * theta_bar;
    for i in 0..6 {
        for j in 0..6 {
            cep[i * 6 + j] -= fac * n[i] * n[j];
        }
    }
    Ok((
        sig,
        cep,
        GpHist {
            ep,
            alpha: hist.alpha + dg,
        },
    ))
}

fn gps(kind: ElemKind, xyz0: &[[f64; 3]]) -> Result<Vec<(Vec<[f64; 3]>, f64)>> {
    match kind {
        ElemKind::Hex8 | ElemKind::Hex8I | ElemKind::Hex8R => {
            let mut p = [[0.0; 3]; 8];
            for i in 0..8 {
                p[i] = xyz0[i];
            }
            let pts = [-G2, G2];
            let mut o = Vec::new();
            for &xi in &pts {
                for &eta in &pts {
                    for &zeta in &pts {
                        let (dndx, det, _) = hex8_dndx(&p, xi, eta, zeta)?;
                        if det <= 0.0 {
                            return err("C3D8 plastic: negative Jakobideterminante.");
                        }
                        o.push((dndx.to_vec(), det));
                    }
                }
            }
            Ok(o)
        }
        ElemKind::Hex20 | ElemKind::Hex20R => {
            let mut o = Vec::new();
            for (xi, eta, zeta, w0) in quadratic::hex_gauss(kind.reduced_int()) {
                let (dndx, det, _) = quadratic::hex20_dndx(xyz0, xi, eta, zeta)?;
                if det <= 0.0 {
                    return err("C3D20 plastic: negative Jakobideterminante.");
                }
                o.push((dndx.to_vec(), w0 * det));
            }
            Ok(o)
        }
        ElemKind::Tet4 => {
            let mut j = [[0.0; 3]; 3];
            for p in 0..3 {
                for q in 0..3 {
                    j[q][p] = xyz0[p + 1][q] - xyz0[0][q];
                }
            }
            let (inv, det) = invert3(j)?;
            if det <= 0.0 {
                return err("C3D4 plastic: negative Jakobideterminante.");
            }
            let vol = det / 6.0;
            let dn = [
                [-1.0, -1.0, -1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ];
            let mut dndx = vec![[0.0; 3]; 4];
            for a in 0..4 {
                for i in 0..3 {
                    dndx[a][i] = inv[0][i] * dn[a][0] + inv[1][i] * dn[a][1] + inv[2][i] * dn[a][2];
                }
            }
            Ok(vec![(dndx, vol)])
        }
        ElemKind::Tet10 | ElemKind::Tet10T => {
            let a = 0.5854101966249685;
            let b = 0.1381966011250105;
            let w0 = 1.0 / 24.0;
            let pts = [[b, b, b], [a, b, b], [b, a, b], [b, b, a]];
            let mut o = Vec::new();
            for p in &pts {
                let (dndx, det, _) = quadratic::tet10_dndx(xyz0, p[0], p[1], p[2])?;
                if det <= 0.0 {
                    return err("C3D10 plastic: negative Jakobideterminante.");
                }
                o.push((dndx.to_vec(), w0 * det));
            }
            Ok(o)
        }
        ElemKind::Wedge6 => {
            let tri = [
                [1.0 / 6.0, 1.0 / 6.0],
                [2.0 / 3.0, 1.0 / 6.0],
                [1.0 / 6.0, 2.0 / 3.0],
            ];
            let wtri = 1.0 / 6.0;
            let mut o = Vec::new();
            for t in &tri {
                for &zeta in &[-G2, G2] {
                    let (dndx, det, _) = crate::extra::wedge_dndx(xyz0, t[0], t[1], zeta)?;
                    if det <= 0.0 {
                        return err("C3D6 plastic: negative Jakobideterminante.");
                    }
                    o.push((dndx.to_vec(), det * wtri));
                }
            }
            Ok(o)
        }
        ElemKind::Wedge15 => {
            let tri = [
                [1.0 / 6.0, 1.0 / 6.0],
                [2.0 / 3.0, 1.0 / 6.0],
                [1.0 / 6.0, 2.0 / 3.0],
            ];
            let wtri = 1.0 / 6.0;
            const G3: [f64; 3] = [-0.7745966692414834, 0.0, 0.7745966692414834];
            const W3: [f64; 3] = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
            let mut o = Vec::new();
            for t in &tri {
                for k in 0..3 {
                    let (dndx, det, _) = crate::extra::wedge15_dndx(xyz0, t[0], t[1], G3[k])?;
                    if det <= 0.0 {
                        return err("C3D15 plastic: negative Jakobideterminante.");
                    }
                    o.push((dndx.to_vec(), det * wtri * W3[k]));
                }
            }
            Ok(o)
        }
        _ => err(format!(
            "*PLASTIC für {} noch nicht implementiert.",
            kind.ccx_name()
        )),
    }
}

pub fn n_gauss(kind: ElemKind) -> usize {
    match kind {
        ElemKind::Hex8 | ElemKind::Hex8I | ElemKind::Hex8R => 8,
        ElemKind::Hex20 => 27,
        ElemKind::Hex20R => 8,
        ElemKind::Tet4 => 1,
        ElemKind::Tet10 | ElemKind::Tet10T => 4,
        ElemKind::Wedge6 => 6,
        ElemKind::Wedge15 => 9,
        _ => 0,
    }
}

/// Assemble tangent and internal force from current displacement, history frozen at step n.
pub fn continuum_plastic(
    kind: ElemKind,
    xyz0: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    curve: &[(f64, f64)],
    hist: &[GpHist],
) -> Result<(PlasticElem, Vec<GpHist>)> {
    let gp = gps(kind, xyz0)?;
    let nn = kind.nnodes();
    let nd = 3 * nn;
    let mut ke = vec![0.0; nd * nd];
    let mut fe = vec![0.0; nd];
    let mut acc_s = [0.0; 6];
    let mut acc_e = [0.0; 6];
    let mut acc_p = 0.0;
    let mut vol = 0.0;
    let mut hist_out = Vec::with_capacity(gp.len());
    let gp = gps(kind, xyz0)?;
    let nn = kind.nnodes();
    let nd = 3 * nn;
    let mut ke = vec![0.0; nd * nd];
    let mut fe = vec![0.0; nd];
    let mut acc_s = [0.0; 6];
    let mut acc_e = [0.0; 6];
    let mut acc_p = 0.0;
    let mut vol = 0.0;
    for (g, (dndx, w)) in gp.iter().enumerate() {
        let mut b = vec![0.0; 6 * nd];
        fill_b3(&mut b, nn, dndx);
        let mut eps = [0.0; 6];
        for i in 0..6 {
            let mut v = 0.0;
            for j in 0..nd {
                v += b[i * nd + j] * ue[j];
            }
            eps[i] = v;
        }
        let h = hist.get(g).copied().unwrap_or_default();
        let (sig, cep, hnew) = j2_return(&eps, &h, e, nu, curve)?;
        gemm_bt_d_b(&mut ke, nd, &b, 6, &cep, *w);
        for j in 0..nd {
            let mut q = 0.0;
            for i in 0..6 {
                q += b[i * nd + j] * sig[i];
            }
            fe[j] += q * *w;
        }
        for i in 0..6 {
            acc_s[i] += sig[i] * *w;
            acc_e[i] += eps[i] * *w;
        }
        acc_p += hnew.alpha * *w;
        vol += *w;
        hist_out.push(hnew);
    }
    let inv = if vol.abs() > 0.0 { 1.0 / vol } else { 0.0 };
    for i in 0..6 {
        acc_s[i] *= inv;
        acc_e[i] *= inv;
    }
    Ok((
        PlasticElem {
            ke,
            fe,
            stress: acc_s,
            strain: acc_e,
            peeq: acc_p * inv,
            ngp: gp.len(),
        },
        hist_out,
    ))
}

/// Total-Lagrange J2: Green–Lagrange strain, PK2 from the small-strain return map,
/// geometric stiffness from S. Small strain recovers `continuum_plastic`.
pub fn continuum_plastic_nl(
    kind: ElemKind,
    xyz0: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    curve: &[(f64, f64)],
    hist: &[GpHist],
) -> Result<(crate::nlgeom::NlElem, Vec<GpHist>)> {
    let gp = gps(kind, xyz0)?;
    let nn = kind.nnodes();
    let nd = 3 * nn;
    let mut ke = vec![0.0; nd * nd];
    let mut fe = vec![0.0; nd];
    let mut acc_s = [0.0; 6];
    let mut acc_e = [0.0; 6];
    let mut acc_p = 0.0;
    let mut vol = 0.0;
    let mut hist_out = Vec::with_capacity(gp.len());
    for (g, (dndx, w)) in gp.iter().enumerate() {
        let f = crate::nlgeom::deformation_gradient(dndx, ue, nn);
        let egl = crate::nlgeom::green_lagrange(&f);
        let eps = [
            egl[0][0],
            egl[1][1],
            egl[2][2],
            2.0 * egl[0][1],
            2.0 * egl[1][2],
            2.0 * egl[2][0],
        ];
        let h = hist.get(g).copied().unwrap_or_default();
        let (s, cep, hnew) = j2_return(&eps, &h, e, nu, curve)?;
        let mut b = vec![0.0; 6 * nd];
        crate::nlgeom::fill_b_nl(&mut b, nn, &f, dndx);
        gemm_bt_d_b(&mut ke, nd, &b, 6, &cep, *w);
        crate::eigen::add_continuum_kg(&mut ke, nn, dndx, &s, *w);
        for j in 0..nd {
            let mut q = 0.0;
            for i in 0..6 {
                q += b[i * nd + j] * s[i];
            }
            fe[j] += q * *w;
        }
        let cauchy = crate::nlgeom::pk2_to_cauchy(&f, &s)?;
        for i in 0..6 {
            acc_s[i] += cauchy[i] * *w;
            acc_e[i] += eps[i] * *w;
        }
        acc_p += hnew.alpha * *w;
        vol += *w;
        hist_out.push(hnew);
    }
    let inv = if vol.abs() > 0.0 { 1.0 / vol } else { 0.0 };
    for i in 0..6 {
        acc_s[i] *= inv;
        acc_e[i] *= inv;
    }
    Ok((
        crate::nlgeom::NlElem {
            ke,
            fe,
            vol,
            cauchy: acc_s,
            gl: acc_e,
            peeq: acc_p * inv,
        },
        hist_out,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn curve() -> Vec<(f64, f64)> {
        vec![(0.0, 210.0), (0.01, 420.0)]
    }

    #[test]
    fn j2_uniaxial_matches_1d() {
        // σ=250, H=21000, peeq=0.00190476, εx=σ/E+peeq=0.00309524
        // εy = −ν σ/E − peeq/2
        let e = 210000.0;
        let nu = 0.3;
        let sig = 250.0;
        let pe = (250.0 - 210.0) / 21000.0;
        let ee = sig / e;
        let eps = [
            ee + pe,
            -nu * ee - 0.5 * pe,
            -nu * ee - 0.5 * pe,
            0.0,
            0.0,
            0.0,
        ];
        let (s, _, h) = j2_return(&eps, &GpHist::default(), e, nu, &curve()).unwrap();
        assert!((s[0] - 250.0).abs() < 0.5, "sxx={}", s[0]);
        assert!(s[1].abs() < 1.0, "syy={}", s[1]);
        assert!(s[2].abs() < 1.0, "szz={}", s[2]);
        assert!((h.alpha - pe).abs() < 1e-6, "peeq={}", h.alpha);
    }

    #[test]
    fn j2_below_yield_is_hooke() {
        let e = 210000.0;
        let nu = 0.3;
        let sig = 200.0;
        let ee = sig / e;
        let eps = [ee, -nu * ee, -nu * ee, 0.0, 0.0, 0.0];
        let (s, _, h) = j2_return(&eps, &GpHist::default(), e, nu, &curve()).unwrap();
        assert!((s[0] - 200.0).abs() < 1e-6, "sxx={}", s[0]);
        assert!(h.alpha.abs() < 1e-16);
    }
}
