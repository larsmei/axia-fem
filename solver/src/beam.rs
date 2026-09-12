//! 3D Timoshenko beam elements (Abaqus B31 / B32).
//!
//! Node order matches Abaqus/CalculiX:
//! - B31: end-1, end-2
//! - B32: end-1, end-2, mid-node
//!
//! Six DOF per node in global axes: u1,u2,u3, ur1,ur2,ur3.
//! Local 1 = tangent t, 2 = n1, 3 = n2 = t × n1.
//! Generalized strains: ε, γ2, γ3, χt, χ2, χ3
//!   γ2 = du_n1/ds − θ_n2
//!   γ3 = du_n2/ds + θ_n1

use crate::error::{err, Result};
use crate::model::{BeamSection, ElemKind};

const G2: f64 = 0.5773502691896257;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn normalize(a: [f64; 3]) -> Result<[f64; 3]> {
    let n = norm(a);
    if n < 1e-18 {
        return err("Balkenrichtung ist degeneriert.");
    }
    Ok(scale(a, 1.0 / n))
}

pub fn orthonormal(t: [f64; 3], n1_hint: [f64; 3]) -> Result<([f64; 3], [f64; 3])> {
    let t = normalize(t)?;
    let mut n1 = add(n1_hint, scale(t, -dot(n1_hint, t)));
    if norm(n1) < 1e-8 {
        let alt = if t[2].abs() < 0.9 {
            [0.0, 0.0, -1.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        n1 = add(alt, scale(t, -dot(alt, t)));
    }
    let n1 = normalize(n1)?;
    let n2 = cross(t, n1);
    Ok((n1, n2))
}

/// Shape N and dN/dξ. ξ ∈ [-1, 1]. Node order: end1, end2 [, mid].
fn shape(kind: ElemKind, xi: f64) -> ([f64; 3], [f64; 3], usize) {
    match kind {
        ElemKind::Beam31 => {
            let n = [(1.0 - xi) * 0.5, (1.0 + xi) * 0.5, 0.0];
            let dn = [-0.5, 0.5, 0.0];
            (n, dn, 2)
        }
        _ => {
            let n = [0.5 * xi * (xi - 1.0), 0.5 * xi * (xi + 1.0), 1.0 - xi * xi];
            let dn = [xi - 0.5, xi + 0.5, -2.0 * xi];
            (n, dn, 3)
        }
    }
}

fn geom_at(xyz: &[[f64; 3]], n: &[f64; 3], dn: &[f64; 3], nn: usize) -> ([f64; 3], [f64; 3], f64) {
    let mut x = [0.0; 3];
    let mut dx = [0.0; 3];
    for a in 0..nn {
        x = add(x, scale(xyz[a], n[a]));
        dx = add(dx, scale(xyz[a], dn[a]));
    }
    let j = norm(dx);
    (x, dx, j)
}

fn gauss(kind: ElemKind) -> &'static [(f64, f64)] {
    match kind {
        ElemKind::Beam31 => &[(0.0, 2.0)],
        _ => &[(-G2, 1.0), (G2, 1.0)],
    }
}

fn dmat(e: f64, nu: f64, sec: &BeamSection) -> Result<[f64; 36]> {
    if e <= 0.0 {
        return err(format!("Ungültiger E-Modul E={e}"));
    }
    if nu >= 0.5 || nu <= -1.0 {
        return err(format!("Ungültige Querkontraktion nu={nu}"));
    }
    let g = e / (2.0 * (1.0 + nu));
    let mut d = [0.0; 36];
    d[0] = e * sec.area;
    d[7] = sec.k11 * g * sec.area;
    d[14] = sec.k22 * g * sec.area;
    d[21] = g * sec.jtor;
    d[28] = e * sec.i11;
    d[35] = e * sec.i22;
    d[28 + 1] = -e * sec.i12; // D[4,5]
    d[35 - 1] = -e * sec.i12; // D[5,4]  index 5*6+4 = 34
    d[34] = -e * sec.i12;
    d[29] = -e * sec.i12;
    Ok(d)
}

/// Fill B (6 × 6nn) mapping *local* generalized strains to *global* nodal dofs.
fn fill_b(
    b: &mut [f64],
    nn: usize,
    nshp: &[f64; 3],
    dnds: &[f64; 3],
    r: &[[f64; 3]; 3], // rows: local basis in global (t, n1, n2)
) {
    // q_local = R^T q_global, R columns = t, n1, n2 so R^T rows = t, n1, n2.
    // r[i][j] = basis_i · e_j  => row i of R^T.
    let nd = 6 * nn;
    b.fill(0.0);
    for a in 0..nn {
        let na = nshp[a];
        let da = dnds[a];
        // Local dofs at node a: [ut, un1, un2, θt, θn1, θn2]
        // Each local component is R^T * (u or θ global)
        // ε  += da * ut
        // γ2 += da * un1 - na * θn2
        // γ3 += da * un2 + na * θn1
        // χt += da * θt
        // χ2 += da * θn1
        // χ3 += da * θn2
        let col0 = a * 6;
        for k in 0..3 {
            // ut contrib from u_g[k]: R^T[0,k] = r[0][k]
            b[0 * nd + col0 + k] += da * r[0][k];
            b[1 * nd + col0 + k] += da * r[1][k];
            b[2 * nd + col0 + k] += da * r[2][k];
            // rotations start at col0+3
            b[1 * nd + col0 + 3 + k] += -na * r[2][k]; // −θn2
            b[2 * nd + col0 + 3 + k] += na * r[1][k]; // +θn1
            b[3 * nd + col0 + 3 + k] += da * r[0][k];
            b[4 * nd + col0 + 3 + k] += da * r[1][k];
            b[5 * nd + col0 + 3 + k] += da * r[2][k];
        }
    }
}

fn gemm_btd_b(ke: &mut [f64], nd: usize, b: &[f64], d: &[f64], w: f64) {
    let nrow = 6usize;
    let mut tmp = [0.0; 6 * 18];
    for i in 0..nrow {
        for j in 0..nd {
            let mut s = 0.0;
            for k in 0..nrow {
                s += d[i * nrow + k] * b[k * nd + j];
            }
            tmp[i * nd + j] = s;
        }
    }
    for i in 0..nd {
        for j in 0..nd {
            let mut s = 0.0;
            for k in 0..nrow {
                s += b[k * nd + i] * tmp[k * nd + j];
            }
            ke[i * nd + j] += w * s;
        }
    }
}

pub fn stiffness(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    e: f64,
    nu: f64,
    sec: &BeamSection,
) -> Result<(Vec<f64>, f64)> {
    let nn = kind.nnodes();
    if xyz.len() < nn {
        return err("Balken hat zu wenige Knotenkoordinaten.");
    }
    let nd = 6 * nn;
    let d = dmat(e, nu, sec)?;
    let mut ke = vec![0.0; nd * nd];
    let mut length = 0.0;
    let mut b = vec![0.0; 6 * nd];
    for &(xi, w) in gauss(kind) {
        let (nshp, dn, _) = shape(kind, xi);
        let (_x, dx, j) = geom_at(xyz, &nshp, &dn, nn);
        if j < 1e-18 {
            return err("Balken-Jacobi ist singulär.");
        }
        let t = scale(dx, 1.0 / j);
        let (n1, n2) = orthonormal(t, sec.n1)?;
        let r = [t, n1, n2];
        let mut dnds = [0.0; 3];
        for a in 0..nn {
            dnds[a] = dn[a] / j;
        }
        fill_b(&mut b, nn, &nshp, &dnds, &r);
        gemm_btd_b(&mut ke, nd, &b, &d, w * j);
        length += w * j;
    }
    Ok((ke, length))
}

fn strains_at(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    ue: &[f64],
    sec: &BeamSection,
    xi: f64,
) -> Result<([f64; 6], [f64; 3], [f64; 3], [f64; 3])> {
    let nn = kind.nnodes();
    let nd = 6 * nn;
    let (nshp, dn, _) = shape(kind, xi);
    let (_x, dx, j) = geom_at(xyz, &nshp, &dn, nn);
    if j < 1e-18 {
        return err("Balken-Jacobi ist singulär.");
    }
    let t = scale(dx, 1.0 / j);
    let (n1, n2) = orthonormal(t, sec.n1)?;
    let r = [t, n1, n2];
    let mut dnds = [0.0; 3];
    for a in 0..nn {
        dnds[a] = dn[a] / j;
    }
    let mut b = vec![0.0; 6 * nd];
    fill_b(&mut b, nn, &nshp, &dnds, &r);
    let mut eps = [0.0; 6];
    for row in 0..6 {
        let mut s = 0.0;
        for jcol in 0..nd.min(ue.len()) {
            s += b[row * nd + jcol] * ue[jcol];
        }
        eps[row] = s;
    }
    Ok((eps, t, n1, n2))
}

/// Nodal Voigt stress from the worst fibre, transformed into global.
pub fn nodal_stress(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    sec: &BeamSection,
) -> Result<Vec<[f64; 6]>> {
    let d = dmat(e, nu, sec)?;
    let nn = kind.nnodes();
    let xis: &[f64] = match kind {
        ElemKind::Beam31 => &[-1.0, 1.0],
        _ => &[-1.0, 1.0, 0.0],
    };
    let mut out = vec![[0.0; 6]; nn];
    for (a, &xi) in xis.iter().enumerate().take(nn) {
        let (eps, t, n1, n2) = strains_at(kind, xyz, ue, sec, xi)?;
        let mut sf = [0.0; 6];
        for i in 0..6 {
            let mut s = 0.0;
            for k in 0..6 {
                s += d[i * 6 + k] * eps[k];
            }
            sf[i] = s;
        }
        // Fibre stress at rectangle/circle extremes in the n1–n2 plane.
        let y = sec.a * 0.5;
        let z = sec.b * 0.5;
        let mut best_s = [0.0; 6];
        let mut best_vm = -1.0;
        for sy in [-1.0, 1.0] {
            for sz in [-1.0, 1.0] {
                let yy = sy * y;
                let zz = sz * z;
                let sx = if sec.area > 0.0 { sf[0] / sec.area } else { 0.0 }
                    - yy * (if sec.i22 > 0.0 { sf[5] / sec.i22 } else { 0.0 })
                    + zz * (if sec.i11 > 0.0 { sf[4] / sec.i11 } else { 0.0 });
                let txy = if sec.area > 0.0 { sf[1] / sec.area } else { 0.0 };
                let txz = if sec.area > 0.0 { sf[2] / sec.area } else { 0.0 };
                let tau_t = if sec.jtor > 0.0 {
                    sf[3] * y.max(z) / sec.jtor
                } else {
                    0.0
                };
                let vm = (sx * sx + 3.0 * (txy * txy + txz * txz + tau_t * tau_t)).sqrt();
                // Global Voigt from σ = sx t⊗t + τ12 (t⊗n1+n1⊗t)/2 + τ13 (t⊗n2+n2⊗t)/2
                let tau12 = txy + tau_t;
                let tau13 = txz;
                let mut sg = [0.0; 6];
                // sxx, syy, szz, sxy, syz, szx
                sg[0] = sx * t[0] * t[0] + 2.0 * tau12 * t[0] * n1[0] + 2.0 * tau13 * t[0] * n2[0];
                sg[1] = sx * t[1] * t[1] + 2.0 * tau12 * t[1] * n1[1] + 2.0 * tau13 * t[1] * n2[1];
                sg[2] = sx * t[2] * t[2] + 2.0 * tau12 * t[2] * n1[2] + 2.0 * tau13 * t[2] * n2[2];
                sg[3] = sx * t[0] * t[1]
                    + tau12 * (t[0] * n1[1] + n1[0] * t[1])
                    + tau13 * (t[0] * n2[1] + n2[0] * t[1]);
                sg[4] = sx * t[1] * t[2]
                    + tau12 * (t[1] * n1[2] + n1[1] * t[2])
                    + tau13 * (t[1] * n2[2] + n2[1] * t[2]);
                sg[5] = sx * t[2] * t[0]
                    + tau12 * (t[2] * n1[0] + n1[2] * t[0])
                    + tau13 * (t[2] * n2[0] + n2[2] * t[0]);
                if vm > best_vm {
                    best_vm = vm;
                    best_s = sg;
                }
            }
        }
        out[a] = best_s;
    }
    Ok(out)
}

fn scatter_trans(
    fe: &mut [f64],
    nn: usize,
    nshp: &[f64; 3],
    wj: f64,
    fx: f64,
    fy: f64,
    fz: f64,
) {
    for a in 0..nn {
        fe[a * 6] += wj * nshp[a] * fx;
        fe[a * 6 + 1] += wj * nshp[a] * fy;
        fe[a * 6 + 2] += wj * nshp[a] * fz;
    }
}

pub fn line_load_global(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    fx: f64,
    fy: f64,
    fz: f64,
) -> Result<Vec<f64>> {
    let nn = kind.nnodes();
    let mut fe = vec![0.0; 6 * nn];
    for &(xi, w) in gauss(kind) {
        let (nshp, dn, _) = shape(kind, xi);
        let (_x, _dx, j) = geom_at(xyz, &nshp, &dn, nn);
        if j < 1e-18 {
            return err("Balken-Jacobi ist singulär.");
        }
        scatter_trans(&mut fe, nn, &nshp, w * j, fx, fy, fz);
    }
    Ok(fe)
}

pub fn body_force(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    sec: &BeamSection,
    bx: f64,
    by: f64,
    bz: f64,
) -> Result<Vec<f64>> {
    line_load_global(kind, xyz, bx * sec.area, by * sec.area, bz * sec.area)
}

pub fn line_load_local(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    sec: &BeamSection,
    axis: usize,
    mag: f64,
) -> Result<Vec<f64>> {
    let nn = kind.nnodes();
    let mut fe = vec![0.0; 6 * nn];
    for &(xi, w) in gauss(kind) {
        let (nshp, dn, _) = shape(kind, xi);
        let (_x, dx, j) = geom_at(xyz, &nshp, &dn, nn);
        if j < 1e-18 {
            return err("Balken-Jacobi ist singulär.");
        }
        let t = scale(dx, 1.0 / j);
        let (n1, n2) = orthonormal(t, sec.n1)?;
        let dir = match axis {
            0 => t,
            1 => n1,
            _ => n2,
        };
        scatter_trans(
            &mut fe,
            nn,
            &nshp,
            w * j,
            mag * dir[0],
            mag * dir[1],
            mag * dir[2],
        );
    }
    Ok(fe)
}
