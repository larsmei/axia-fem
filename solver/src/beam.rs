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

/// Nastran CBAR: Hermitian bending with Timoshenko shear, exact on a straight member.
/// Local DOF per node: u_t, u_n1, u_n2, θ_t, θ_n1, θ_n2.
pub fn cbar_stiffness(
    xyz: &[[f64; 3]],
    e: f64,
    nu: f64,
    sec: &BeamSection,
) -> Result<(Vec<f64>, f64)> {
    if xyz.len() < 2 {
        return err("CBAR braucht zwei Knoten.");
    }
    if e <= 0.0 {
        return err(format!("Ungültiger E-Modul E={e}"));
    }
    let p0 = [
        xyz[0][0] + sec.off_a[0],
        xyz[0][1] + sec.off_a[1],
        xyz[0][2] + sec.off_a[2],
    ];
    let p1 = [
        xyz[1][0] + sec.off_b[0],
        xyz[1][1] + sec.off_b[1],
        xyz[1][2] + sec.off_b[2],
    ];
    let dx = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let len = norm(dx);
    if len < 1e-18 {
        return err("CBAR hat die Länge null.");
    }
    let t = scale(dx, 1.0 / len);
    let (n1, n2) = orthonormal(t, sec.n1)?;
    let g = e / (2.0 * (1.0 + nu).max(1e-6));
    let mut kl = [0.0; 144];
    let ea = e * sec.area / len;
    kl[0 * 12 + 0] = ea;
    kl[0 * 12 + 6] = -ea;
    kl[6 * 12 + 0] = -ea;
    kl[6 * 12 + 6] = ea;
    let gj = g * sec.jtor / len;
    kl[3 * 12 + 3] = gj;
    kl[3 * 12 + 9] = -gj;
    kl[9 * 12 + 3] = -gj;
    kl[9 * 12 + 9] = gj;
    // Deflection along n1, rotation about n2 (θ_euler = θ_n2), inertia I22, shear k11.
    add_timoshenko(&mut kl, e * sec.i22, sec.k11 * g * sec.area, len, [1, 5, 7, 11], [1.0, 1.0, 1.0, 1.0]);
    // Deflection along n2, θ_euler = −θ_n1, inertia I11, shear k22.
    add_timoshenko(&mut kl, e * sec.i11, sec.k22 * g * sec.area, len, [2, 4, 8, 10], [1.0, -1.0, 1.0, -1.0]);
    condense_releases(&mut kl, sec.rel_a, sec.rel_b);
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        r[i][0] = t[i];
        r[i][1] = n1[i];
        r[i][2] = n2[i];
    }
    let mut ke = vec![0.0; 144];
    // u_g = R u_l, K_g = T K_l T^T
    let mut tl = [0.0; 144];
    for node in 0..2 {
        for blk in 0..2 {
            let o = node * 6 + blk * 3;
            for i in 0..3 {
                for j in 0..3 {
                    tl[(o + i) * 12 + (o + j)] = r[i][j];
                }
            }
        }
    }
    let mut tmp = [0.0; 144];
    for i in 0..12 {
        for j in 0..12 {
            let mut s = 0.0;
            for k in 0..12 {
                s += tl[i * 12 + k] * kl[k * 12 + j];
            }
            tmp[i * 12 + j] = s;
        }
    }
    for i in 0..12 {
        for j in 0..12 {
            let mut s = 0.0;
            for k in 0..12 {
                s += tmp[i * 12 + k] * tl[j * 12 + k];
            }
            ke[i * 12 + j] = s;
        }
    }
    apply_end_offsets(&mut ke, sec.off_a, sec.off_b);
    Ok((ke, len))
}

fn add_timoshenko(kl: &mut [f64], ei: f64, ks: f64, len: f64, idx: [usize; 4], sign: [f64; 4]) {
    if ei.abs() < 1e-30 {
        return;
    }
    let phi = if ks.abs() < 1e-18 {
        0.0
    } else {
        (12.0 * ei / (ks * len * len)).clamp(0.0, 1.0e6)
    };
    let c = ei / ((1.0 + phi) * len * len * len);
    let l = len;
    let kel = [
        [12.0 * c, 6.0 * l * c, -12.0 * c, 6.0 * l * c],
        [6.0 * l * c, (4.0 + phi) * l * l * c, -6.0 * l * c, (2.0 - phi) * l * l * c],
        [-12.0 * c, -6.0 * l * c, 12.0 * c, -6.0 * l * c],
        [6.0 * l * c, (2.0 - phi) * l * l * c, -6.0 * l * c, (4.0 + phi) * l * l * c],
    ];
    for a in 0..4 {
        for b in 0..4 {
            let v = sign[a] * kel[a][b] * sign[b];
            kl[idx[a] * 12 + idx[b]] += v;
        }
    }
}

fn condense_releases(kl: &mut [f64], rel_a: u8, rel_b: u8) {
    let mut free = Vec::new();
    for i in 0..6 {
        if rel_a & (1 << i) != 0 {
            free.push(i);
        }
        if rel_b & (1 << i) != 0 {
            free.push(6 + i);
        }
    }
    for &s in &free {
        let kss = kl[s * 12 + s];
        if kss.abs() < 1e-18 {
            continue;
        }
        for i in 0..12 {
            if i == s {
                continue;
            }
            let kis = kl[i * 12 + s];
            if kis.abs() == 0.0 {
                continue;
            }
            for j in 0..12 {
                if j == s {
                    continue;
                }
                kl[i * 12 + j] -= kis * kl[s * 12 + j] / kss;
            }
        }
        for i in 0..12 {
            kl[i * 12 + s] = 0.0;
            kl[s * 12 + i] = 0.0;
        }
    }
}

fn apply_end_offsets(ke: &mut [f64], wa: [f64; 3], wb: [f64; 3]) {
    let za = wa[0].abs() + wa[1].abs() + wa[2].abs();
    let zb = wb[0].abs() + wb[1].abs() + wb[2].abs();
    if za + zb < 1e-18 {
        return;
    }
    // u_beam = u_grid + θ_grid × w, θ_beam = θ_grid. w points from the grid to the beam end.
    let mut t = [0.0; 144];
    for n in 0..2 {
        let w = if n == 0 { wa } else { wb };
        let o = n * 6;
        for i in 0..6 {
            t[(o + i) * 12 + (o + i)] = 1.0;
        }
        let s = [
            [0.0, w[2], -w[1]],
            [-w[2], 0.0, w[0]],
            [w[1], -w[0], 0.0],
        ];
        for i in 0..3 {
            for j in 0..3 {
                t[(o + i) * 12 + (o + 3 + j)] = s[i][j];
            }
        }
    }
    let mut tmp = [0.0; 144];
    for i in 0..12 {
        for j in 0..12 {
            let mut s = 0.0;
            for k in 0..12 {
                s += ke[i * 12 + k] * t[k * 12 + j];
            }
            tmp[i * 12 + j] = s;
        }
    }
    for i in 0..12 {
        for j in 0..12 {
            let mut s = 0.0;
            for k in 0..12 {
                s += t[k * 12 + i] * tmp[k * 12 + j];
            }
            ke[i * 12 + j] = s;
        }
    }
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

/// Co-rotational NLGEOM: linear Timoshenko in the current chord frame.
pub fn stiffness_nl(
    kind: ElemKind,
    xyz0: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    sec: &BeamSection,
) -> Result<(Vec<f64>, Vec<f64>, [f64; 6])> {
    let nn = kind.nnodes();
    let nd = 6 * nn;
    let mut xyz = vec![[0.0; 3]; nn];
    for a in 0..nn {
        xyz[a] = [
            xyz0[a][0] + ue.get(6 * a).copied().unwrap_or(0.0),
            xyz0[a][1] + ue.get(6 * a + 1).copied().unwrap_or(0.0),
            xyz0[a][2] + ue.get(6 * a + 2).copied().unwrap_or(0.0),
        ];
    }
    let i1 = if nn == 2 { 1 } else { 1 }; // B32: nodes 0,1 are ends (order: end1,end2,mid)
    let i1 = if nn == 3 { 1 } else { i1 };
    let t0 = normalize([
        xyz0[i1][0] - xyz0[0][0],
        xyz0[i1][1] - xyz0[0][1],
        xyz0[i1][2] - xyz0[0][2],
    ])?;
    let t = normalize([
        xyz[i1][0] - xyz[0][0],
        xyz[i1][1] - xyz[0][1],
        xyz[i1][2] - xyz[0][2],
    ])?;
    let l0 = dist3(xyz0[0], xyz0[i1]).max(1e-18);
    let l = dist3(xyz[0], xyz[i1]).max(1e-18);
    let (n1, n2) = orthonormal(t, sec.n1)?;
    let r = rotation_from_t(t0, t);
    let th_chord = rotvec_from_r(r);
    let mut udef = vec![0.0; nd];
    for a in 0..nn {
        let mut thg = [
            ue.get(6 * a + 3).copied().unwrap_or(0.0),
            ue.get(6 * a + 4).copied().unwrap_or(0.0),
            ue.get(6 * a + 5).copied().unwrap_or(0.0),
        ];
        thg = [
            thg[0] - th_chord[0],
            thg[1] - th_chord[1],
            thg[2] - th_chord[2],
        ];
        udef[6 * a + 3] = t[0] * thg[0] + t[1] * thg[1] + t[2] * thg[2];
        udef[6 * a + 4] = n1[0] * thg[0] + n1[1] * thg[1] + n1[2] * thg[2];
        udef[6 * a + 5] = n2[0] * thg[0] + n2[1] * thg[1] + n2[2] * thg[2];
    }
    udef[6 * i1] = l - l0;
    if nn == 3 {
        udef[12] = 0.5 * (l - l0);
    }
    let (ke_l, _) = stiffness(kind, xyz0, e, nu, sec)?;
    let mut fe_l = vec![0.0; nd];
    for i in 0..nd {
        let mut s = 0.0;
        for j in 0..nd {
            s += ke_l[i * nd + j] * udef[j];
        }
        fe_l[i] = s;
    }
    let mut ke = vec![0.0; nd * nd];
    let mut fe = vec![0.0; nd];
    rotate_6(nn, &ke_l, &fe_l, t, n1, n2, &mut ke, &mut fe);
    let nforce = e * sec.area * (l - l0) / l0;
    let geom = nforce / l;
    for a in [0usize, i1] {
        for b in [0usize, i1] {
            let sg = if a == b { geom } else { -geom };
            for i in 0..3 {
                for j in 0..3 {
                    let pr = if i == j { 1.0 } else { 0.0 } - t[i] * t[j];
                    ke[(6 * a + i) * nd + (6 * b + j)] += sg * pr;
                }
            }
        }
    }
    let sig = if sec.area > 0.0 { nforce / sec.area } else { 0.0 };
    let cauchy = [
        sig * t[0] * t[0],
        sig * t[1] * t[1],
        sig * t[2] * t[2],
        sig * t[0] * t[1],
        sig * t[1] * t[2],
        sig * t[2] * t[0],
    ];
    Ok((ke, fe, cauchy))
}

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

fn rotation_from_t(t0: [f64; 3], t: [f64; 3]) -> [[f64; 3]; 3] {
    let c = dot(t0, t).clamp(-1.0, 1.0);
    let mut axis = cross(t0, t);
    let s = norm(axis);
    if s < 1e-14 {
        if c >= 0.0 {
            return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        }
        axis = if t0[2].abs() < 0.9 {
            normalize(cross(t0, [0.0, 0.0, 1.0])).unwrap_or([1.0, 0.0, 0.0])
        } else {
            normalize(cross(t0, [0.0, 1.0, 0.0])).unwrap_or([1.0, 0.0, 0.0])
        };
        let e = axis;
        return [
            [2.0 * e[0] * e[0] - 1.0, 2.0 * e[0] * e[1], 2.0 * e[0] * e[2]],
            [2.0 * e[1] * e[0], 2.0 * e[1] * e[1] - 1.0, 2.0 * e[1] * e[2]],
            [2.0 * e[2] * e[0], 2.0 * e[2] * e[1], 2.0 * e[2] * e[2] - 1.0],
        ];
    }
    axis = scale(axis, 1.0 / s);
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        r[i][i] += c;
        r[i][(i + 1) % 3] += -axis[(i + 2) % 3] * s;
        r[i][(i + 2) % 3] += axis[(i + 1) % 3] * s;
        for j in 0..3 {
            r[i][j] += (1.0 - c) * axis[i] * axis[j];
        }
    }
    r
}

fn rotvec_from_r(r: [[f64; 3]; 3]) -> [f64; 3] {
    let c = ((r[0][0] + r[1][1] + r[2][2] - 1.0) * 0.5).clamp(-1.0, 1.0);
    let ang = c.acos();
    if ang.abs() < 1e-14 {
        return [0.0, 0.0, 0.0];
    }
    let s = ang.sin();
    if s.abs() < 1e-14 {
        return [0.0, 0.0, 0.0];
    }
    [
        ang * (r[2][1] - r[1][2]) / (2.0 * s),
        ang * (r[0][2] - r[2][0]) / (2.0 * s),
        ang * (r[1][0] - r[0][1]) / (2.0 * s),
    ]
}

fn rotate_6(
    nn: usize,
    ke_l: &[f64],
    fe_l: &[f64],
    t: [f64; 3],
    n1: [f64; 3],
    n2: [f64; 3],
    ke: &mut [f64],
    fe: &mut [f64],
) {
    let nd = 6 * nn;
    let q = [
        [t[0], n1[0], n2[0]],
        [t[1], n1[1], n2[1]],
        [t[2], n1[2], n2[2]],
    ];
    let mut tfull = vec![0.0; nd * nd];
    for a in 0..nn {
        for blk in 0..2 {
            let o = 6 * a + 3 * blk;
            for i in 0..3 {
                for j in 0..3 {
                    tfull[(o + i) * nd + (o + j)] = q[i][j];
                }
            }
        }
    }
    let mut tmp = vec![0.0; nd * nd];
    for i in 0..nd {
        for j in 0..nd {
            let mut s = 0.0;
            for k in 0..nd {
                s += tfull[i * nd + k] * ke_l[k * nd + j];
            }
            tmp[i * nd + j] = s;
        }
    }
    for i in 0..nd {
        for j in 0..nd {
            let mut s = 0.0;
            for k in 0..nd {
                s += tmp[i * nd + k] * tfull[j * nd + k];
            }
            ke[i * nd + j] = s;
        }
        let mut s = 0.0;
        for k in 0..nd {
            s += tfull[i * nd + k] * fe_l[k];
        }
        fe[i] = s;
    }
}

