//! Total-Lagrange St. Venant–Kirchhoff for continuum (NLGEOM).
//! Small strain recovers linear elasticity.

use crate::elem::{d_iso_3d, gemm_bt_d_b, hex8_dndx, invert3, G2};
use crate::eigen::add_continuum_kg;
use crate::error::{err, Result};
use crate::model::ElemKind;
use crate::quadratic;

pub struct NlElem {
    pub ke: Vec<f64>,
    pub fe: Vec<f64>,
    #[allow(dead_code)]
    pub vol: f64,
    pub cauchy: [f64; 6],
    pub gl: [f64; 6],
    pub peeq: f64,
}

fn voigt_s(e: &[[f64; 3]; 3], lam: f64, mu: f64) -> [f64; 6] {
    let tr = e[0][0] + e[1][1] + e[2][2];
    [
        lam * tr + 2.0 * mu * e[0][0],
        lam * tr + 2.0 * mu * e[1][1],
        lam * tr + 2.0 * mu * e[2][2],
        2.0 * mu * e[0][1],
        2.0 * mu * e[1][2],
        2.0 * mu * e[2][0],
    ]
}

pub(crate) fn green_lagrange(f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = f[0][i] * f[0][j] + f[1][i] * f[1][j] + f[2][i] * f[2][j];
        }
    }
    let mut e = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            e[i][j] = 0.5 * (c[i][j] - if i == j { 1.0 } else { 0.0 });
        }
    }
    e
}

pub(crate) fn deformation_gradient(dndx: &[[f64; 3]], ue: &[f64], nn: usize) -> [[f64; 3]; 3] {
    let mut f = [[0.0; 3]; 3];
    for i in 0..3 {
        f[i][i] = 1.0;
    }
    for a in 0..nn {
        for i in 0..3 {
            let u = ue[3 * a + i];
            for j in 0..3 {
                f[i][j] += u * dndx[a][j];
            }
        }
    }
    f
}

fn det3(a: &[[f64; 3]; 3]) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

pub(crate) fn pk2_to_cauchy(f: &[[f64; 3]; 3], s: &[f64; 6]) -> Result<[f64; 6]> {
    let sm = [
        [s[0], s[3], s[5]],
        [s[3], s[1], s[4]],
        [s[5], s[4], s[2]],
    ];
    let mut fs = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            fs[i][j] = f[i][0] * sm[0][j] + f[i][1] * sm[1][j] + f[i][2] * sm[2][j];
        }
    }
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = fs[i][0] * f[j][0] + fs[i][1] * f[j][1] + fs[i][2] * f[j][2];
        }
    }
    let j = det3(f);
    if j <= 1e-18 {
        return err("NLGEOM: det(F) ≤ 0 (Element inversion).");
    }
    let invj = 1.0 / j;
    Ok([
        c[0][0] * invj,
        c[1][1] * invj,
        c[2][2] * invj,
        c[0][1] * invj,
        c[1][2] * invj,
        c[2][0] * invj,
    ])
}

pub(crate) fn fill_b_nl(b: &mut [f64], nn: usize, f: &[[f64; 3]; 3], dndx: &[[f64; 3]]) {
    // Voigt E_eng: [Exx, Eyy, Ezz, 2Exy, 2Eyz, 2Ezx]
    let n = 3 * nn;
    for a in 0..nn {
        for i in 0..3 {
            let c = 3 * a + i;
            b[0 * n + c] = f[i][0] * dndx[a][0];
            b[1 * n + c] = f[i][1] * dndx[a][1];
            b[2 * n + c] = f[i][2] * dndx[a][2];
            b[3 * n + c] = f[i][0] * dndx[a][1] + f[i][1] * dndx[a][0];
            b[4 * n + c] = f[i][1] * dndx[a][2] + f[i][2] * dndx[a][1];
            b[5 * n + c] = f[i][2] * dndx[a][0] + f[i][0] * dndx[a][2];
        }
    }
}

fn lam_mu(e: f64, nu: f64) -> Result<(f64, f64)> {
    if e <= 0.0 || nu >= 0.5 || nu <= -1.0 {
        return err(format!("Ungültiges Material E={e} nu={nu}"));
    }
    let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    Ok((lam, mu))
}

/// Tangent, internal force, reference volume, volume-averaged Cauchy and Green–Lagrange.
pub fn continuum_nl(
    kind: ElemKind,
    xyz0: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
) -> Result<NlElem> {
    match kind {
        ElemKind::Hex8 | ElemKind::Hex8I | ElemKind::Hex8R => hex8_nl(xyz0, ue, e, nu),
        ElemKind::Hex20 | ElemKind::Hex20R => hex20_nl(xyz0, ue, e, nu, kind.reduced_int()),
        ElemKind::Tet4 => tet4_nl(xyz0, ue, e, nu),
        ElemKind::Tet10 | ElemKind::Tet10T => tet10_nl(xyz0, ue, e, nu),
        ElemKind::Wedge6 => wedge6_nl(xyz0, ue, e, nu),
        ElemKind::Wedge15 => wedge15_nl(xyz0, ue, e, nu),
        _ => err(format!(
            "NLGEOM für {} noch nicht implementiert.",
            kind.ccx_name()
        )),
    }
}

struct GpAcc {
    ke: Vec<f64>,
    fe: Vec<f64>,
    vol: f64,
    cauchy: [f64; 6],
    gl: [f64; 6],
    peeq: f64,
}

impl GpAcc {
    fn new(nd: usize) -> Self {
        Self {
            ke: vec![0.0; nd * nd],
            fe: vec![0.0; nd],
            vol: 0.0,
            cauchy: [0.0; 6],
            gl: [0.0; 6],
            peeq: 0.0,
        }
    }
    fn finish(self) -> NlElem {
        let inv = if self.vol.abs() > 0.0 {
            1.0 / self.vol
        } else {
            0.0
        };
        let mut cauchy = self.cauchy;
        let mut gl = self.gl;
        for i in 0..6 {
            cauchy[i] *= inv;
            gl[i] *= inv;
        }
        NlElem {
            ke: self.ke,
            fe: self.fe,
            vol: self.vol,
            cauchy,
            gl,
            peeq: self.peeq * inv,
        }
    }
}

fn assemble_gp(
    acc: &mut GpAcc,
    nn: usize,
    dndx: &[[f64; 3]],
    ue: &[f64],
    cmat: &[f64],
    lam: f64,
    mu: f64,
    w: f64,
) -> Result<()> {
    let nd = 3 * nn;
    let f = deformation_gradient(dndx, ue, nn);
    if det3(&f) <= 1e-18 {
        return err("NLGEOM: det(F) ≤ 0 (Element inversion).");
    }
    let egl = green_lagrange(&f);
    let s = voigt_s(&egl, lam, mu);
    let mut b = vec![0.0; 6 * nd];
    fill_b_nl(&mut b, nn, &f, dndx);
    gemm_bt_d_b(&mut acc.ke, nd, &b, 6, cmat, w);
    add_continuum_kg(&mut acc.ke, nn, dndx, &s, w);
    for j in 0..nd {
        let mut q = 0.0;
        for i in 0..6 {
            q += b[i * nd + j] * s[i];
        }
        acc.fe[j] += q * w;
    }
    let cauchy = pk2_to_cauchy(&f, &s)?;
    let glv = [
        egl[0][0], egl[1][1], egl[2][2], 2.0 * egl[0][1], 2.0 * egl[1][2], 2.0 * egl[2][0],
    ];
    for i in 0..6 {
        acc.cauchy[i] += cauchy[i] * w;
        acc.gl[i] += glv[i] * w;
    }
    acc.vol += w;
    Ok(())
}

fn hex8_nl(xyz0: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<NlElem> {
    let mut p = [[0.0; 3]; 8];
    for i in 0..8 {
        p[i] = xyz0[i];
    }
    let (lam, mu) = lam_mu(e, nu)?;
    let cmat = d_iso_3d(e, nu)?;
    let mut acc = GpAcc::new(24);
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            for &zeta in &pts {
                let (dndx, det, _) = hex8_dndx(&p, xi, eta, zeta)?;
                if det <= 0.0 {
                    return err("NLGEOM C3D8: negative Jakobideterminante.");
                }
                assemble_gp(&mut acc, 8, &dndx, ue, &cmat, lam, mu, det)?;
            }
        }
    }
    Ok(acc.finish())
}

fn hex20_nl(
    xyz0: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    reduced: bool,
) -> Result<NlElem> {
    let (lam, mu) = lam_mu(e, nu)?;
    let cmat = d_iso_3d(e, nu)?;
    let mut acc = GpAcc::new(60);
    for (xi, eta, zeta, w0) in quadratic::hex_gauss(reduced) {
        let (dndx, det, _) = quadratic::hex20_dndx(xyz0, xi, eta, zeta)?;
        if det <= 0.0 {
            return err("NLGEOM C3D20: negative Jakobideterminante.");
        }
        assemble_gp(&mut acc, 20, &dndx, ue, &cmat, lam, mu, w0 * det)?;
    }
    Ok(acc.finish())
}

fn tet4_nl(xyz0: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<NlElem> {
    let mut j = [[0.0; 3]; 3];
    for p in 0..3 {
        for q in 0..3 {
            j[q][p] = xyz0[p + 1][q] - xyz0[0][q];
        }
    }
    let (inv, det) = invert3(j)?;
    if det <= 0.0 {
        return err("NLGEOM C3D4: negative Jakobideterminante.");
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
    let (lam, mu) = lam_mu(e, nu)?;
    let cmat = d_iso_3d(e, nu)?;
    let mut acc = GpAcc::new(12);
    assemble_gp(&mut acc, 4, &dndx, ue, &cmat, lam, mu, vol)?;
    Ok(acc.finish())
}

fn tet10_nl(xyz0: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<NlElem> {
    let (lam, mu) = lam_mu(e, nu)?;
    let cmat = d_iso_3d(e, nu)?;
    let mut acc = GpAcc::new(30);
    let a = 0.5854101966249685;
    let b = 0.1381966011250105;
    let w0 = 1.0 / 24.0;
    let pts = [[b, b, b], [a, b, b], [b, a, b], [b, b, a]];
    for p in &pts {
        let (dndx, det, _) = quadratic::tet10_dndx(xyz0, p[0], p[1], p[2])?;
        if det <= 0.0 {
            return err("NLGEOM C3D10: negative Jakobideterminante.");
        }
        assemble_gp(&mut acc, 10, &dndx, ue, &cmat, lam, mu, w0 * det)?;
    }
    Ok(acc.finish())
}

fn wedge6_nl(xyz0: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<NlElem> {
    let (lam, mu) = lam_mu(e, nu)?;
    let cmat = d_iso_3d(e, nu)?;
    let mut acc = GpAcc::new(18);
    let tri = [
        [1.0 / 6.0, 1.0 / 6.0],
        [2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0],
    ];
    let wtri = 1.0 / 6.0;
    for t in &tri {
        for &zeta in &[-G2, G2] {
            let (dndx, det, _) = crate::extra::wedge_dndx(xyz0, t[0], t[1], zeta)?;
            if det <= 0.0 {
                return err("NLGEOM C3D6: negative Jakobideterminante.");
            }
            assemble_gp(&mut acc, 6, &dndx, ue, &cmat, lam, mu, det * wtri)?;
        }
    }
    Ok(acc.finish())
}

fn wedge15_nl(xyz0: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<NlElem> {
    let (lam, mu) = lam_mu(e, nu)?;
    let cmat = d_iso_3d(e, nu)?;
    let mut acc = GpAcc::new(45);
    let tri = [
        [1.0 / 6.0, 1.0 / 6.0],
        [2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0],
    ];
    let wtri = 1.0 / 6.0;
    const G3: [f64; 3] = [-0.7745966692414834, 0.0, 0.7745966692414834];
    const W3: [f64; 3] = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
    for t in &tri {
        for k in 0..3 {
            let (dndx, det, _) = crate::extra::wedge15_dndx(xyz0, t[0], t[1], G3[k])?;
            if det <= 0.0 {
                return err("NLGEOM C3D15: negative Jakobideterminante.");
            }
            assemble_gp(
                &mut acc,
                15,
                &dndx,
                ue,
                &cmat,
                lam,
                mu,
                det * wtri * W3[k],
            )?;
        }
    }
    Ok(acc.finish())
}

pub fn is_nl_continuum(kind: ElemKind) -> bool {
    kind.is_continuum3d()
}
