//! Total-Lagrange St. Venant–Kirchhoff for continuum (NLGEOM).
//! Small strain recovers linear elasticity.

use crate::eigen::add_continuum_kg;
use crate::elem::{
    d_iso_3d, d_plane_strain, d_plane_stress, gemm_bt_d_b, hex8_dndx, invert2, invert3, G2,
};
use crate::error::{err, Result};
use crate::model::{ElemKind, HyperKind, Material};
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
    let sm = [[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]];
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
    mat: &Material,
    thickness: f64,
) -> Result<NlElem> {
    if kind.is_axisym() {
        return cax_nl(kind, xyz0, ue, mat);
    }
    if is_nl_plane(kind) {
        return plane_nl(kind, xyz0, ue, mat, thickness);
    }
    match kind {
        ElemKind::Hex8 | ElemKind::Hex8I | ElemKind::Hex8R => hex8_nl(xyz0, ue, mat),
        ElemKind::Hex20 | ElemKind::Hex20R => hex20_nl(xyz0, ue, mat, kind.reduced_int()),
        ElemKind::Tet4 => tet4_nl(xyz0, ue, mat),
        ElemKind::Tet10 | ElemKind::Tet10T => tet10_nl(xyz0, ue, mat),
        ElemKind::Wedge6 => wedge6_nl(xyz0, ue, mat),
        ElemKind::Wedge15 => wedge15_nl(xyz0, ue, mat),
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
    mat: &Material,
    w: f64,
) -> Result<()> {
    let nd = 3 * nn;
    let f = deformation_gradient(dndx, ue, nn);
    if det3(&f) <= 1e-18 {
        return err("NLGEOM: det(F) ≤ 0 (Element inversion).");
    }
    let (s, cmat) = pk2_and_c(&f, mat)?;
    let mut b = vec![0.0; 6 * nd];
    fill_b_nl(&mut b, nn, &f, dndx);
    gemm_bt_d_b(&mut acc.ke, nd, &b, 6, &cmat, w);
    add_continuum_kg(&mut acc.ke, nn, dndx, &s, w);
    for j in 0..nd {
        let mut q = 0.0;
        for i in 0..6 {
            q += b[i * nd + j] * s[i];
        }
        acc.fe[j] += q * w;
    }
    let cauchy = pk2_to_cauchy(&f, &s)?;
    let egl = green_lagrange(&f);
    let glv = [
        egl[0][0],
        egl[1][1],
        egl[2][2],
        2.0 * egl[0][1],
        2.0 * egl[1][2],
        2.0 * egl[2][0],
    ];
    for i in 0..6 {
        acc.cauchy[i] += cauchy[i] * w;
        acc.gl[i] += glv[i] * w;
    }
    acc.vol += w;
    Ok(())
}

fn hex8_nl(xyz0: &[[f64; 3]], ue: &[f64], mat: &Material) -> Result<NlElem> {
    let mut p = [[0.0; 3]; 8];
    for i in 0..8 {
        p[i] = xyz0[i];
    }
    let mut acc = GpAcc::new(24);
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            for &zeta in &pts {
                let (dndx, det, _) = hex8_dndx(&p, xi, eta, zeta)?;
                if det <= 0.0 {
                    return err("NLGEOM C3D8: negative Jakobideterminante.");
                }
                assemble_gp(&mut acc, 8, &dndx, ue, mat, det)?;
            }
        }
    }
    Ok(acc.finish())
}

fn hex20_nl(xyz0: &[[f64; 3]], ue: &[f64], mat: &Material, reduced: bool) -> Result<NlElem> {
    let mut acc = GpAcc::new(60);
    for (xi, eta, zeta, w0) in quadratic::hex_gauss(reduced) {
        let (dndx, det, _) = quadratic::hex20_dndx(xyz0, xi, eta, zeta)?;
        if det <= 0.0 {
            return err("NLGEOM C3D20: negative Jakobideterminante.");
        }
        assemble_gp(&mut acc, 20, &dndx, ue, mat, w0 * det)?;
    }
    Ok(acc.finish())
}

/// Cauchy and Green–Lagrange at the 20 C3D20 nodes (Gauss + Lagrange, like ccx).
pub fn hex20_nodal_nl(
    xyz0: &[[f64; 3]],
    ue: &[f64],
    mat: &Material,
    reduced: bool,
) -> Result<(Vec<[f64; 6]>, Vec<[f64; 6]>)> {
    let gps = quadratic::hex_gauss(reduced);
    let mut g_c = Vec::with_capacity(gps.len());
    let mut g_e = Vec::with_capacity(gps.len());
    for &(xi, eta, zeta, _) in &gps {
        let (dndx, det, _) = quadratic::hex20_dndx(xyz0, xi, eta, zeta)?;
        if det <= 0.0 {
            return err("NLGEOM C3D20: negative Jakobideterminante.");
        }
        let f = deformation_gradient(&dndx, ue, 20);
        if det3(&f) <= 1e-18 {
            return err("NLGEOM: det(F) ≤ 0 (Element inversion).");
        }
        let (s, _) = pk2_and_c(&f, mat)?;
        g_c.push(pk2_to_cauchy(&f, &s)?);
        let egl = green_lagrange(&f);
        g_e.push([
            egl[0][0],
            egl[1][1],
            egl[2][2],
            2.0 * egl[0][1],
            2.0 * egl[1][2],
            2.0 * egl[2][0],
        ]);
    }
    Ok((
        quadratic::hex20_extrapolate(&g_c, reduced),
        quadratic::hex20_extrapolate(&g_e, reduced),
    ))
}

fn tet4_nl(xyz0: &[[f64; 3]], ue: &[f64], mat: &Material) -> Result<NlElem> {
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
    let mut acc = GpAcc::new(12);
    assemble_gp(&mut acc, 4, &dndx, ue, mat, vol)?;
    Ok(acc.finish())
}

fn tet10_nl(xyz0: &[[f64; 3]], ue: &[f64], mat: &Material) -> Result<NlElem> {
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
        assemble_gp(&mut acc, 10, &dndx, ue, mat, w0 * det)?;
    }
    Ok(acc.finish())
}

fn wedge6_nl(xyz0: &[[f64; 3]], ue: &[f64], mat: &Material) -> Result<NlElem> {
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
            assemble_gp(&mut acc, 6, &dndx, ue, mat, det * wtri)?;
        }
    }
    Ok(acc.finish())
}

fn wedge15_nl(xyz0: &[[f64; 3]], ue: &[f64], mat: &Material) -> Result<NlElem> {
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
            assemble_gp(&mut acc, 15, &dndx, ue, mat, det * wtri * W3[k])?;
        }
    }
    Ok(acc.finish())
}

pub fn is_nl_plane(kind: ElemKind) -> bool {
    matches!(
        kind,
        ElemKind::Quad4Ps
            | ElemKind::Quad4Pe
            | ElemKind::Quad8Ps
            | ElemKind::Quad8Pe
            | ElemKind::Quad8RPs
            | ElemKind::Quad8RPe
            | ElemKind::Tri3Ps
            | ElemKind::Tri3Pe
            | ElemKind::Tri6Ps
            | ElemKind::Tri6Pe
    )
}

pub fn is_nl_continuum(kind: ElemKind) -> bool {
    kind.is_continuum3d() || is_nl_plane(kind) || kind.is_axisym()
}

fn pk2_and_c(f: &[[f64; 3]; 3], mat: &Material) -> Result<([f64; 6], [f64; 36])> {
    if mat.is_hyper() {
        return hyper_pk2_and_c(f, mat);
    }
    let egl = green_lagrange(f);
    let (lam, mu) = lam_mu(mat.e, mat.nu)?;
    let s = voigt_s(&egl, lam, mu);
    let cmat = d_iso_3d(mat.e, mat.nu)?;
    Ok((s, cmat))
}

fn c_from_f(f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = f[0][i] * f[0][j] + f[1][i] * f[1][j] + f[2][i] * f[2][j];
        }
    }
    c
}

fn inv3_det(a: &[[f64; 3]; 3]) -> Result<([[f64; 3]; 3], f64)> {
    invert3(*a)
}

fn pk2_hyper(c: &[[f64; 3]; 3], mat: &Material) -> Result<[f64; 6]> {
    let (cinv, detc) = inv3_det(c)?;
    if detc <= 1e-30 {
        return err("NLGEOM hyperelastisch: det(C) ≤ 0.");
    }
    let j = detc.sqrt();
    match mat.hyper {
        HyperKind::NeoHooke => {
            let c10 = mat.h[0];
            let d1 = mat.h[1];
            let mu = 2.0 * c10;
            let kappa = if d1.abs() > 1e-18 {
                2.0 / d1
            } else {
                1.0e4 * mu.abs().max(1.0)
            };
            let i1 = c[0][0] + c[1][1] + c[2][2];
            let jm23 = j.powf(-2.0 / 3.0);
            let ibar = i1 * jm23;
            let vol = kappa * (j - 1.0) * j;
            let iso = mu * jm23;
            let mut s = [[0.0; 3]; 3];
            for i in 0..3 {
                for k in 0..3 {
                    let ident = if i == k { 1.0 } else { 0.0 };
                    s[i][k] = iso * (ident - ibar / 3.0 * cinv[i][k]) + vol * cinv[i][k];
                }
            }
            Ok([s[0][0], s[1][1], s[2][2], s[0][1], s[1][2], s[2][0]])
        }
        HyperKind::Ogden | HyperKind::Hyperfoam => ogden_pk2(c, j, &cinv, mat),
        HyperKind::None => err("kein Hyperelast"),
    }
}

fn ogden_pk2(c: &[[f64; 3]; 3], j: f64, cinv: &[[f64; 3]; 3], mat: &Material) -> Result<[f64; 6]> {
    let (evals, evecs) = eigen3_sym(c)?;
    let mut lam = [0.0; 3];
    for i in 0..3 {
        lam[i] = evals[i].max(1e-16).sqrt();
    }
    let nterm = mat.h_n.max(1) as usize;
    let mut t_prin = [0.0; 3]; // Kirchhoff τ_i / (later PK2)
    for t in 0..nterm.min(2) {
        let (mu, alpha, beta_or_d) = if mat.hyper == HyperKind::Hyperfoam {
            (mat.h[3 * t], mat.h[3 * t + 1], mat.h[3 * t + 2])
        } else {
            // Ogden packed μ1,α1,D1 [,μ2,α2,D2] — D only on first term
            if t == 0 {
                (mat.h[0], mat.h[1], mat.h[2])
            } else {
                (mat.h[3], mat.h[4], mat.h[5])
            }
        };
        if mu.abs() < 1e-18 || alpha.abs() < 1e-18 {
            continue;
        }
        for i in 0..3 {
            t_prin[i] += mu * (lam[i].powf(alpha) - 1.0);
        }
        if mat.hyper == HyperKind::Hyperfoam {
            let beta = beta_or_d;
            if beta.abs() > 1e-18 {
                let p = mu / alpha * (-alpha * beta) * j.powf(-alpha * beta);
                for i in 0..3 {
                    t_prin[i] += p;
                }
            }
        }
    }
    if mat.hyper == HyperKind::Ogden {
        let d1 = mat.h[2];
        let kappa = if d1.abs() > 1e-18 {
            2.0 / d1
        } else {
            1.0e4 * mat.h[0].abs().max(1.0)
        };
        let p = kappa * (j - 1.0) * j;
        for i in 0..3 {
            t_prin[i] += p;
        }
        let _ = cinv;
    }
    // τ_i = t_prin (Kirchhoff-like), S = F^{-1} τ F^{-T} = sum (τ_i / λ_i²) n⊗n
    let mut s = [[0.0; 3]; 3];
    for i in 0..3 {
        let coef = t_prin[i] / evals[i].max(1e-16);
        for a in 0..3 {
            for b in 0..3 {
                s[a][b] += coef * evecs[a][i] * evecs[b][i];
            }
        }
    }
    Ok([s[0][0], s[1][1], s[2][2], s[0][1], s[1][2], s[2][0]])
}

fn eigen3_sym(a: &[[f64; 3]; 3]) -> Result<([f64; 3], [[f64; 3]; 3])> {
    // Jacobi
    let mut s = *a;
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..32 {
        let mut p = 0.0;
        let mut pi = 0;
        let mut pj = 1;
        for i in 0..3 {
            for j in (i + 1)..3 {
                if s[i][j].abs() > p {
                    p = s[i][j].abs();
                    pi = i;
                    pj = j;
                }
            }
        }
        if p < 1e-14 {
            break;
        }
        let app = s[pi][pi];
        let aqq = s[pj][pj];
        let apq = s[pi][pj];
        let tau = 0.5 * (aqq - app) / apq;
        let t = if tau.abs() > 1e12 {
            0.5 / tau
        } else {
            tau.signum() / (tau.abs() + (1.0 + tau * tau).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let ss = t * c;
        s[pi][pi] = app - t * apq;
        s[pj][pj] = aqq + t * apq;
        s[pi][pj] = 0.0;
        s[pj][pi] = 0.0;
        for k in 0..3 {
            if k != pi && k != pj {
                let aik = s[pi][k];
                let ajk = s[pj][k];
                s[pi][k] = c * aik - ss * ajk;
                s[k][pi] = s[pi][k];
                s[pj][k] = ss * aik + c * ajk;
                s[k][pj] = s[pj][k];
            }
            let vip = v[k][pi];
            let viq = v[k][pj];
            v[k][pi] = c * vip - ss * viq;
            v[k][pj] = ss * vip + c * viq;
        }
    }
    Ok(([s[0][0], s[1][1], s[2][2]], v))
}

fn hyper_pk2_and_c(f: &[[f64; 3]; 3], mat: &Material) -> Result<([f64; 6], [f64; 36])> {
    let c0 = c_from_f(f);
    let s0 = pk2_hyper(&c0, mat)?;
    let egl = green_lagrange(f);
    let e0 = [
        egl[0][0],
        egl[1][1],
        egl[2][2],
        2.0 * egl[0][1],
        2.0 * egl[1][2],
        2.0 * egl[2][0],
    ];
    let eps = 1e-7;
    let mut d = [0.0; 36];
    for k in 0..6 {
        let mut ep = e0;
        ep[k] += eps;
        let mut cp = [
            [1.0 + 2.0 * ep[0], ep[3], ep[5]],
            [ep[3], 1.0 + 2.0 * ep[1], ep[4]],
            [ep[5], ep[4], 1.0 + 2.0 * ep[2]],
        ];
        // symmetrize
        cp[1][0] = cp[0][1];
        cp[2][0] = cp[0][2];
        cp[2][1] = cp[1][2];
        let sp = pk2_hyper(&cp, mat).unwrap_or(s0);
        for i in 0..6 {
            d[i * 6 + k] = (sp[i] - s0[i]) / eps;
        }
    }
    Ok((s0, d))
}

fn plane_nl(
    kind: ElemKind,
    xyz0: &[[f64; 3]],
    ue: &[f64],
    mat: &Material,
    thickness: f64,
) -> Result<NlElem> {
    let nn = kind.nnodes();
    let th = thickness.abs().max(1e-12);
    let plane_strain = kind.is_plane_strain();
    let reduced = kind.reduced_int();
    let mut acc = GpAcc::new(3 * nn);
    for (xi, eta, w0) in plane_gauss(nn, reduced) {
        let (dndx2, det, _) = plane_dndx(kind, xyz0, xi, eta)?;
        if det <= 0.0 {
            return err(format!(
                "NLGEOM {}: negative Jakobideterminante.",
                kind.ccx_name()
            ));
        }
        let mut dndx3 = vec![[0.0; 3]; nn];
        for a in 0..nn {
            dndx3[a] = [dndx2[a][0], dndx2[a][1], 0.0];
        }
        if plane_strain || mat.is_hyper() {
            assemble_gp(&mut acc, nn, &dndx3, ue, mat, th * w0 * det)?;
        } else {
            assemble_gp_plane_stress(&mut acc, nn, &dndx2, ue, mat, th * w0 * det)?;
        }
    }
    Ok(acc.finish())
}

fn assemble_gp_plane_stress(
    acc: &mut GpAcc,
    nn: usize,
    dndx2: &[[f64; 2]],
    ue: &[f64],
    mat: &Material,
    w: f64,
) -> Result<()> {
    let nd = 3 * nn;
    let mut f2 = [[1.0, 0.0], [0.0, 1.0]];
    for a in 0..nn {
        let ux = ue[3 * a];
        let uy = ue[3 * a + 1];
        f2[0][0] += ux * dndx2[a][0];
        f2[0][1] += ux * dndx2[a][1];
        f2[1][0] += uy * dndx2[a][0];
        f2[1][1] += uy * dndx2[a][1];
    }
    let detf = f2[0][0] * f2[1][1] - f2[0][1] * f2[1][0];
    if detf <= 1e-18 {
        return err("NLGEOM plane stress: det(F) ≤ 0.");
    }
    let egl = [
        [
            0.5 * (f2[0][0] * f2[0][0] + f2[1][0] * f2[1][0] - 1.0),
            0.5 * (f2[0][0] * f2[0][1] + f2[1][0] * f2[1][1]),
        ],
        [0.0, 0.5 * (f2[0][1] * f2[0][1] + f2[1][1] * f2[1][1] - 1.0)],
    ];
    let e_eng = [egl[0][0], egl[1][1], 2.0 * egl[0][1]];
    let d = d_plane_stress(mat.e, mat.nu)?;
    let mut s2 = [0.0; 3];
    for i in 0..3 {
        for k in 0..3 {
            s2[i] += d[i * 3 + k] * e_eng[k];
        }
    }
    let mut b = vec![0.0; 3 * nd];
    for a in 0..nn {
        for i in 0..2 {
            let c = 3 * a + i;
            b[0 * nd + c] = f2[i][0] * dndx2[a][0];
            b[1 * nd + c] = f2[i][1] * dndx2[a][1];
            b[2 * nd + c] = f2[i][0] * dndx2[a][1] + f2[i][1] * dndx2[a][0];
        }
    }
    gemm_bt_d_b(&mut acc.ke, nd, &b, 3, &d, w);
    // geometric stiffness in-plane
    let sxx = s2[0];
    let syy = s2[1];
    let sxy = s2[2];
    for a in 0..nn {
        for bnode in 0..nn {
            let ga = dndx2[a];
            let gb = dndx2[bnode];
            let gtg = ga[0] * (sxx * gb[0] + sxy * gb[1]) + ga[1] * (sxy * gb[0] + syy * gb[1]);
            let v = gtg * w;
            for dir in 0..2 {
                acc.ke[(3 * a + dir) * nd + (3 * bnode + dir)] += v;
            }
        }
    }
    for j in 0..nd {
        let mut q = 0.0;
        for i in 0..3 {
            q += b[i * nd + j] * s2[i];
        }
        acc.fe[j] += q * w;
    }
    let invj = 1.0 / detf;
    acc.cauchy[0] += (f2[0][0] * (s2[0] * f2[0][0] + s2[2] * f2[0][1])
        + f2[0][1] * (s2[2] * f2[0][0] + s2[1] * f2[0][1]))
        * invj
        * w;
    acc.cauchy[1] += (f2[1][0] * (s2[0] * f2[1][0] + s2[2] * f2[1][1])
        + f2[1][1] * (s2[2] * f2[1][0] + s2[1] * f2[1][1]))
        * invj
        * w;
    acc.cauchy[3] += (f2[0][0] * (s2[0] * f2[1][0] + s2[2] * f2[1][1])
        + f2[0][1] * (s2[2] * f2[1][0] + s2[1] * f2[1][1]))
        * invj
        * w;
    acc.gl[0] += e_eng[0] * w;
    acc.gl[1] += e_eng[1] * w;
    acc.gl[3] += e_eng[2] * w;
    acc.vol += w;
    let _ = d_plane_strain;
    Ok(())
}

fn plane_gauss(nn: usize, reduced: bool) -> Vec<(f64, f64, f64)> {
    if nn == 3 {
        return vec![(1.0 / 3.0, 1.0 / 3.0, 0.5)];
    }
    if nn == 6 {
        let a = 1.0 / 6.0;
        let b = 2.0 / 3.0;
        return vec![(a, a, 1.0 / 6.0), (b, a, 1.0 / 6.0), (a, b, 1.0 / 6.0)];
    }
    if reduced && nn <= 4 {
        return vec![(0.0, 0.0, 4.0)];
    }
    if reduced || nn == 4 {
        let g = G2;
        let mut o = Vec::new();
        for &xi in &[-g, g] {
            for &eta in &[-g, g] {
                o.push((xi, eta, 1.0));
            }
        }
        return o;
    }
    const G3: [f64; 3] = [-0.7745966692414834, 0.0, 0.7745966692414834];
    const W3: [f64; 3] = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
    let mut o = Vec::new();
    for i in 0..3 {
        for j in 0..3 {
            o.push((G3[i], G3[j], W3[i] * W3[j]));
        }
    }
    o
}

fn plane_dndx(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    xi: f64,
    eta: f64,
) -> Result<(Vec<[f64; 2]>, f64, Vec<f64>)> {
    let nn = kind.nnodes();
    let mut xy = vec![[0.0; 2]; nn];
    for a in 0..nn {
        xy[a] = [xyz[a][0], xyz[a][1]];
    }
    match nn {
        4 => {
            let mut p = [[0.0; 2]; 4];
            p.copy_from_slice(&xy[..4]);
            let (dndx, det, n) = quad4_dndx_local(&p, xi, eta)?;
            Ok((dndx.to_vec(), det, n.to_vec()))
        }
        8 => {
            let (dndx, det, n) = quadratic::quad8_dndx(&xy, xi, eta)?;
            Ok((dndx.to_vec(), det, n.to_vec()))
        }
        3 => {
            let (dndx, two_a) = tri3_dndx_local(&xy);
            if two_a <= 0.0 {
                return err("NLGEOM CPS3/CPE3: Fläche ≤ 0.");
            }
            Ok((dndx.to_vec(), two_a, vec![1.0 / 3.0; 3]))
        }
        6 => {
            let (dndx, det, n) = quadratic::tri6_dndx(&xy, xi, eta)?;
            Ok((dndx.to_vec(), det, n.to_vec()))
        }
        _ => err("NLGEOM plane: unbekannte Knotenzahl"),
    }
}

fn quad4_dndx_local(
    xy: &[[f64; 2]; 4],
    xi: f64,
    eta: f64,
) -> Result<([[f64; 2]; 4], f64, [f64; 4])> {
    let mut n = [0.0; 4];
    let mut dn = [[0.0; 2]; 4];
    const Q: [[f64; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
    for i in 0..4 {
        n[i] = 0.25 * (1.0 + Q[i][0] * xi) * (1.0 + Q[i][1] * eta);
        dn[i][0] = 0.25 * Q[i][0] * (1.0 + Q[i][1] * eta);
        dn[i][1] = 0.25 * Q[i][1] * (1.0 + Q[i][0] * xi);
    }
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

fn tri3_dndx_local(xy: &[[f64; 2]]) -> ([[f64; 2]; 3], f64) {
    let x1 = xy[0][0];
    let y1 = xy[0][1];
    let x2 = xy[1][0];
    let y2 = xy[1][1];
    let x3 = xy[2][0];
    let y3 = xy[2][1];
    let two_a = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    let inv = if two_a.abs() > 1e-18 {
        1.0 / two_a
    } else {
        0.0
    };
    let dndx = [
        [(y2 - y3) * inv, (x3 - x2) * inv],
        [(y3 - y1) * inv, (x1 - x3) * inv],
        [(y1 - y2) * inv, (x2 - x1) * inv],
    ];
    (dndx, two_a)
}

fn cax_nl(kind: ElemKind, xyz0: &[[f64; 3]], ue: &[f64], mat: &Material) -> Result<NlElem> {
    let nn = kind.nnodes();
    let reduced = kind.reduced_int();
    let mut acc = GpAcc::new(3 * nn);
    for (xi, eta, w0) in plane_gauss(nn, reduced) {
        let (dndx2, det, nshp) = plane_dndx(kind, xyz0, xi, eta)?;
        if det <= 0.0 {
            return err(format!(
                "NLGEOM {}: negative Jakobideterminante.",
                kind.ccx_name()
            ));
        }
        let mut r0 = 0.0;
        for a in 0..nn {
            r0 += nshp[a] * xyz0[a][0];
        }
        let mut ur = 0.0;
        for a in 0..nn {
            ur += nshp[a] * ue[3 * a];
        }
        let r = (r0 + ur).abs().max(1e-12);
        let mut dndx3 = vec![[0.0; 3]; nn];
        for a in 0..nn {
            dndx3[a] = [dndx2[a][0], dndx2[a][1], 0.0];
        }
        // Hoop stretch Fθθ = r/R is injected by adding a dummy dN/dθ via a 3rd axis
        // using an equivalent ∇N_θ = N/R on the radial displacement.
        // We build F explicitly and a specialized B.
        let w = 2.0 * std::f64::consts::PI * r0.abs().max(1e-12) * w0 * det;
        assemble_gp_cax(&mut acc, nn, &dndx2, &nshp, r0, ue, mat, w, r)?;
        let _ = dndx3;
    }
    Ok(acc.finish())
}

fn assemble_gp_cax(
    acc: &mut GpAcc,
    nn: usize,
    dndx2: &[[f64; 2]],
    nshp: &[f64],
    r0: f64,
    ue: &[f64],
    mat: &Material,
    w: f64,
    r: f64,
) -> Result<()> {
    let nd = 3 * nn;
    let mut f = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, r / r0.abs().max(1e-12)],
    ];
    for a in 0..nn {
        let ur = ue[3 * a];
        let uz = ue[3 * a + 1];
        f[0][0] += ur * dndx2[a][0];
        f[0][1] += ur * dndx2[a][1];
        f[1][0] += uz * dndx2[a][0];
        f[1][1] += uz * dndx2[a][1];
    }
    if det3(&f) <= 1e-18 {
        return err("NLGEOM CAX: det(F) ≤ 0.");
    }
    let (s, cmat) = pk2_and_c(&f, mat)?;
    let mut dndx3 = vec![[0.0; 3]; nn];
    let invr = 1.0 / r0.abs().max(1e-12);
    for a in 0..nn {
        dndx3[a] = [dndx2[a][0], dndx2[a][1], nshp[a] * invr];
    }
    let mut b = vec![0.0; 6 * nd];
    fill_b_nl(&mut b, nn, &f, &dndx3);
    gemm_bt_d_b(&mut acc.ke, nd, &b, 6, &cmat, w);
    add_continuum_kg(&mut acc.ke, nn, &dndx3, &s, w);
    for j in 0..nd {
        let mut q = 0.0;
        for i in 0..6 {
            q += b[i * nd + j] * s[i];
        }
        acc.fe[j] += q * w;
    }
    let cauchy = pk2_to_cauchy(&f, &s)?;
    let egl = green_lagrange(&f);
    let glv = [
        egl[0][0],
        egl[1][1],
        egl[2][2],
        2.0 * egl[0][1],
        2.0 * egl[1][2],
        2.0 * egl[2][0],
    ];
    for i in 0..6 {
        acc.cauchy[i] += cauchy[i] * w;
        acc.gl[i] += glv[i] * w;
    }
    acc.vol += w;
    let _ = r;
    Ok(())
}
