//! Reissner–Mindlin shells, Abaqus-compatible:
//! S4 / S4R (MITC4), S3 (DKT + CST), S8 / S8R (SRI), S6.
//!
//! 6 DOF per node in global axes: u1,u2,u3, ur1,ur2,ur3.
//! Local: membrane plane-stress + bending + transverse shear + drilling.

use crate::elem::{d_plane_stress, fill_b2, gemm_bt_d_b, invert2, G2};
use crate::error::{err, Result};
use crate::model::ElemKind;
use crate::quadratic::{quad8_shape, tri6_shape};

const G3: [f64; 3] = [-0.7745966692414834, 0.0, 0.7745966692414834];
const W3: [f64; 3] = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
const K_SHEAR: f64 = 5.0 / 6.0;
const DRILL: f64 = 1.0e-4;

const QUAD_XI: [[f64; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];

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
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn normalize(a: [f64; 3]) -> Result<[f64; 3]> {
    let n = norm(a);
    if n < 1e-18 {
        return err("Schalennormale ist degeneriert.");
    }
    Ok(scale(a, 1.0 / n))
}

/// Orthonormal triad: e1, e2 in-plane, e3 = normal (right-handed, CCW nodes).
pub fn local_frame(xyz: &[[f64; 3]], nn: usize) -> Result<([f64; 3], [f64; 3], [f64; 3])> {
    let e3 = if nn >= 4 {
        let gxi = scale(add(sub(xyz[1], xyz[0]), sub(xyz[2], xyz[3])), 0.25);
        let geta = scale(add(sub(xyz[3], xyz[0]), sub(xyz[2], xyz[1])), 0.25);
        let n = cross(gxi, geta);
        if norm(n) > 1e-14 {
            normalize(n)?
        } else {
            normalize(cross(sub(xyz[1], xyz[0]), sub(xyz[3], xyz[0])))?
        }
    } else {
        normalize(cross(sub(xyz[1], xyz[0]), sub(xyz[2], xyz[0])))?
    };
    let mut e1 = sub(xyz[1], xyz[0]);
    e1 = sub(e1, scale(e3, dot(e1, e3)));
    if norm(e1) < 1e-12 {
        let alt = if e3[2].abs() < 0.9 {
            [0.0, 0.0, 1.0]
        } else {
            [1.0, 0.0, 0.0]
        };
        e1 = sub(alt, scale(e3, dot(alt, e3)));
    }
    let e1 = normalize(e1)?;
    let e2 = cross(e3, e1);
    Ok((e1, e2, e3))
}

fn project_xy(xyz: &[[f64; 3]], e1: [f64; 3], e2: [f64; 3], nn: usize) -> Vec<[f64; 2]> {
    let o = xyz[0];
    (0..nn)
        .map(|i| {
            let d = sub(xyz[i], o);
            [dot(d, e1), dot(d, e2)]
        })
        .collect()
}

fn t6(e1: [f64; 3], e2: [f64; 3], e3: [f64; 3]) -> [[f64; 6]; 6] {
    let mut t = [[0.0; 6]; 6];
    for i in 0..3 {
        t[i][0] = e1[i];
        t[i][1] = e2[i];
        t[i][2] = e3[i];
        t[i + 3][3] = e1[i];
        t[i + 3][4] = e2[i];
        t[i + 3][5] = e3[i];
    }
    t
}

fn rotate_ke(ke: &mut [f64], nn: usize, e1: [f64; 3], e2: [f64; 3], e3: [f64; 3]) {
    let nd = 6 * nn;
    let t = t6(e1, e2, e3);
    let mut kg = vec![0.0; nd * nd];
    for a in 0..nn {
        for b in 0..nn {
            let mut kl = [[0.0; 6]; 6];
            for i in 0..6 {
                for j in 0..6 {
                    kl[i][j] = ke[(6 * a + i) * nd + (6 * b + j)];
                }
            }
            let mut klt = [[0.0; 6]; 6];
            for i in 0..6 {
                for j in 0..6 {
                    for k in 0..6 {
                        klt[i][j] += kl[i][k] * t[j][k];
                    }
                }
            }
            for i in 0..6 {
                for j in 0..6 {
                    let mut s = 0.0;
                    for k in 0..6 {
                        s += t[i][k] * klt[k][j];
                    }
                    kg[(6 * a + i) * nd + (6 * b + j)] = s;
                }
            }
        }
    }
    ke.copy_from_slice(&kg);
}

fn rotate_vec_to_local(ug: &[f64], nn: usize, e1: [f64; 3], e2: [f64; 3], e3: [f64; 3]) -> Vec<f64> {
    let t = t6(e1, e2, e3);
    let mut ul = vec![0.0; 6 * nn];
    for a in 0..nn {
        for i in 0..6 {
            let mut s = 0.0;
            for k in 0..6 {
                s += t[k][i] * ug[6 * a + k];
            }
            ul[6 * a + i] = s;
        }
    }
    ul
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

fn jac_xy(xy: &[[f64; 2]], dn: &[[f64; 2]], nn: usize) -> Result<([[f64; 2]; 2], f64, Vec<[f64; 2]>)> {
    let mut j = [[0.0; 2]; 2];
    for a in 0..nn {
        j[0][0] += dn[a][0] * xy[a][0];
        j[0][1] += dn[a][1] * xy[a][0];
        j[1][0] += dn[a][0] * xy[a][1];
        j[1][1] += dn[a][1] * xy[a][1];
    }
    let (inv, det) = invert2(j)?;
    let mut dndx = vec![[0.0; 2]; nn];
    for a in 0..nn {
        dndx[a][0] = inv[0][0] * dn[a][0] + inv[1][0] * dn[a][1];
        dndx[a][1] = inv[0][1] * dn[a][0] + inv[1][1] * dn[a][1];
    }
    Ok((j, det, dndx))
}

fn add_membrane(ke: &mut [f64], nd: usize, nn: usize, dndx: &[[f64; 2]], dm: &[f64], w: f64) {
    let n = 2 * nn;
    let mut b = vec![0.0; 3 * n];
    fill_b2(&mut b, nn, dndx);
    let mut ke2 = vec![0.0; n * n];
    gemm_bt_d_b(&mut ke2, n, &b, 3, dm, w);
    for a in 0..nn {
        for b in 0..nn {
            for i in 0..2 {
                for j in 0..2 {
                    ke[(6 * a + i) * nd + (6 * b + j)] += ke2[(2 * a + i) * n + (2 * b + j)];
                }
            }
        }
    }
}

fn add_bending(ke: &mut [f64], nd: usize, nn: usize, dndx: &[[f64; 2]], db: &[f64], w: f64) {
    // κ = [θy,x,  -θx,y,  θy,y - θx,x]
    let n = 2 * nn; // θx, θy packed
    let mut b = vec![0.0; 3 * n];
    for i in 0..nn {
        let c = 2 * i;
        b[0 * n + c + 1] = dndx[i][0]; // κx ← θy,x
        b[1 * n + c] = -dndx[i][1]; // κy ← -θx,y
        b[2 * n + c] = -dndx[i][0]; // κxy ← -θx,x
        b[2 * n + c + 1] = dndx[i][1]; // κxy ← θy,y
    }
    let mut ke2 = vec![0.0; n * n];
    gemm_bt_d_b(&mut ke2, n, &b, 3, db, w);
    for a in 0..nn {
        for b in 0..nn {
            for i in 0..2 {
                for j in 0..2 {
                    ke[(6 * a + 3 + i) * nd + (6 * b + 3 + j)] += ke2[(2 * a + i) * n + (2 * b + j)];
                }
            }
        }
    }
}

fn add_shear(ke: &mut [f64], nd: usize, bg: &[f64], ds: &[f64], w: f64, nn: usize) {
    // bg is 2 x (6*nn), γ = [γxz, γyz]
    let n = 6 * nn;
    let mut ke2 = vec![0.0; n * n];
    gemm_bt_d_b(&mut ke2, n, bg, 2, ds, w);
    for i in 0..n {
        for j in 0..n {
            ke[i * nd + j] += ke2[i * n + j];
        }
    }
}

fn ds_mat(e: f64, nu: f64, h: f64) -> [f64; 4] {
    let g = e / (2.0 * (1.0 + nu)) * K_SHEAR * h;
    [g, 0.0, 0.0, g]
}

fn add_drill(ke: &mut [f64], nd: usize, nn: usize, e: f64, h: f64, area: f64) {
    let k = DRILL * e * h * area / nn as f64;
    for a in 0..nn {
        let i = 6 * a + 5;
        ke[i * nd + i] += k;
    }
}

/// MITC4 covariant shear B (2 × 24) at (ξ,η).
fn mitc4_bgamma(xy: &[[f64; 2]], xi: f64, eta: f64) -> Result<(Vec<f64>, f64)> {
    let nn = 4usize;
    let n = 24usize;
    let mut bcov = vec![0.0; 2 * n];
    // γ_ξζ from tying (0, ±1)
    for (eta_t, coeff) in [(1.0, 0.5 * (1.0 + eta)), (-1.0, 0.5 * (1.0 - eta))] {
        let (nshp, dn) = quad4_shape(0.0, eta_t);
        let mut xxi = 0.0;
        let mut yxi = 0.0;
        for a in 0..nn {
            xxi += dn[a][0] * xy[a][0];
            yxi += dn[a][0] * xy[a][1];
        }
        for a in 0..nn {
            // γ_ξζ = w,ξ + x,ξ θy − y,ξ θx
            bcov[0 * n + 6 * a + 2] += coeff * dn[a][0];
            bcov[0 * n + 6 * a + 3] += coeff * (-nshp[a] * yxi);
            bcov[0 * n + 6 * a + 4] += coeff * (nshp[a] * xxi);
        }
    }
    for (xi_t, coeff) in [(1.0, 0.5 * (1.0 + xi)), (-1.0, 0.5 * (1.0 - xi))] {
        let (nshp, dn) = quad4_shape(xi_t, 0.0);
        let mut xet = 0.0;
        let mut yet = 0.0;
        for a in 0..nn {
            xet += dn[a][1] * xy[a][0];
            yet += dn[a][1] * xy[a][1];
        }
        for a in 0..nn {
            bcov[1 * n + 6 * a + 2] += coeff * dn[a][1];
            bcov[1 * n + 6 * a + 3] += coeff * (-nshp[a] * yet);
            bcov[1 * n + 6 * a + 4] += coeff * (nshp[a] * xet);
        }
    }
    let (_, dn) = quad4_shape(xi, eta);
    let (j, det, _) = jac_xy(xy, &dn, nn)?;
    // [γxz; γyz] = J^{-T} [γ_ξζ; γ_ηζ]
    let (inv, _) = invert2(j)?;
    // J^{-T}[i][k] = inv[k][i]
    let mut bg = vec![0.0; 2 * n];
    for col in 0..n {
        let gxi = bcov[col];
        let get = bcov[n + col];
        bg[col] = inv[0][0] * gxi + inv[1][0] * get;
        bg[n + col] = inv[0][1] * gxi + inv[1][1] * get;
    }
    Ok((bg, det))
}

fn cartesian_shear_b(nshp: &[f64], dndx: &[[f64; 2]], nn: usize) -> Vec<f64> {
    // γxz = w,x + θy ; γyz = w,y − θx
    let n = 6 * nn;
    let mut b = vec![0.0; 2 * n];
    for a in 0..nn {
        b[0 * n + 6 * a + 2] = dndx[a][0];
        b[0 * n + 6 * a + 4] = nshp[a];
        b[1 * n + 6 * a + 2] = dndx[a][1];
        b[1 * n + 6 * a + 3] = -nshp[a];
    }
    b
}

fn s4_local(xy: &[[f64; 2]], e: f64, nu: f64, h: f64) -> Result<(Vec<f64>, f64)> {
    let nn = 4usize;
    let nd = 24usize;
    let mut ke = vec![0.0; nd * nd];
    let dm0 = d_plane_stress(e, nu)?;
    let mut dm = [0.0; 9];
    let mut db = [0.0; 9];
    for i in 0..9 {
        dm[i] = dm0[i] * h;
        db[i] = dm0[i] * h * h * h / 12.0;
    }
    let ds = ds_mat(e, nu, h);
    let mut area = 0.0;
    for &xi in &[-G2, G2] {
        for &eta in &[-G2, G2] {
            let (nshp, dn) = quad4_shape(xi, eta);
            let (_, det, dndx) = jac_xy(xy, &dn, nn)?;
            if det <= 0.0 {
                return err("S4: negative Jakobideterminante.");
            }
            add_membrane(&mut ke, nd, nn, &dndx, &dm, det);
            add_bending(&mut ke, nd, nn, &dndx, &db, det);
            let (bg, _) = mitc4_bgamma(xy, xi, eta)?;
            add_shear(&mut ke, nd, &bg, &ds, det, nn);
            let _ = nshp;
            area += det;
        }
    }
    add_drill(&mut ke, nd, nn, e, h, area);
    Ok((ke, area))
}

fn s8_local(xy: &[[f64; 2]], e: f64, nu: f64, h: f64, reduced_shear: bool) -> Result<(Vec<f64>, f64)> {
    let nn = 8usize;
    let nd = 48usize;
    let mut ke = vec![0.0; nd * nd];
    let dm0 = d_plane_stress(e, nu)?;
    let mut dm = [0.0; 9];
    let mut db = [0.0; 9];
    for i in 0..9 {
        dm[i] = dm0[i] * h;
        db[i] = dm0[i] * h * h * h / 12.0;
    }
    let ds = ds_mat(e, nu, h);
    let mut area = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let xi = G3[i];
            let eta = G3[j];
            let w = W3[i] * W3[j];
            let (nshp, dn) = quad8_shape(xi, eta);
            let mut dna = [[0.0; 2]; 8];
            dna.copy_from_slice(&dn);
            let (_, det, dndx) = jac_xy(xy, &dna, nn)?;
            if det <= 0.0 {
                return err("S8: negative Jakobideterminante.");
            }
            add_membrane(&mut ke, nd, nn, &dndx, &dm, w * det);
            add_bending(&mut ke, nd, nn, &dndx, &db, w * det);
            if !reduced_shear {
                let bg = cartesian_shear_b(&nshp, &dndx, nn);
                add_shear(&mut ke, nd, &bg, &ds, w * det, nn);
            }
            area += w * det;
        }
    }
    if reduced_shear {
        for &xi in &[-G2, G2] {
            for &eta in &[-G2, G2] {
                let (nshp, dn) = quad8_shape(xi, eta);
                let mut dna = [[0.0; 2]; 8];
                dna.copy_from_slice(&dn);
                let (_, det, dndx) = jac_xy(xy, &dna, nn)?;
                let bg = cartesian_shear_b(&nshp, &dndx, nn);
                add_shear(&mut ke, nd, &bg, &ds, det, nn);
            }
        }
    }
    add_drill(&mut ke, nd, nn, e, h, area);
    Ok((ke, area))
}

fn tri3_shape() -> ([f64; 3], [[f64; 2]; 3]) {
    // parent: L2=ξ, L3=η, L1=1-ξ-η ; used only for pressure/body. dN in parent.
    (
        [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]],
    )
}

/// DKT bending 9×9 (w, θx, θy) plus CST membrane, 1-pt shear.
fn s3_local(xy: &[[f64; 2]; 3], e: f64, nu: f64, h: f64) -> Result<(Vec<f64>, f64)> {
    let x1 = xy[0][0];
    let y1 = xy[0][1];
    let x2 = xy[1][0];
    let y2 = xy[1][1];
    let x3 = xy[2][0];
    let y3 = xy[2][1];
    let two_a = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
    if two_a <= 0.0 {
        return err("S3: nicht-positive Fläche.");
    }
    let area = 0.5 * two_a;
    let nd = 18usize;
    let mut ke = vec![0.0; nd * nd];

    // CST membrane
    let mut dndx = [[0.0; 2]; 3];
    dndx[0] = [(y2 - y3) / two_a, (x3 - x2) / two_a];
    dndx[1] = [(y3 - y1) / two_a, (x1 - x3) / two_a];
    dndx[2] = [(y1 - y2) / two_a, (x2 - x1) / two_a];
    let dm0 = d_plane_stress(e, nu)?;
    let mut dm = [0.0; 9];
    let mut db = [0.0; 9];
    for i in 0..9 {
        dm[i] = dm0[i] * h;
        db[i] = dm0[i] * h * h * h / 12.0;
    }
    add_membrane(&mut ke, nd, 3, &dndx, &dm, area);

    // DKT bending, 3 Hammer points
    let x23 = x2 - x3;
    let y23 = y2 - y3;
    let x31 = x3 - x1;
    let y31 = y3 - y1;
    let x12 = x1 - x2;
    let y12 = y1 - y2;
    let l23 = x23 * x23 + y23 * y23;
    let l31 = x31 * x31 + y31 * y31;
    let l12 = x12 * x12 + y12 * y12;
    let p4 = -6.0 * y23 / l23;
    let p5 = -6.0 * y31 / l31;
    let p6 = -6.0 * y12 / l12;
    let q4 = 3.0 * y23 * y23 / l23;
    let q5 = 3.0 * y31 * y31 / l31;
    let q6 = 3.0 * y12 * y12 / l12;
    let r4 = 3.0 * x23 * y23 / l23;
    let r5 = 3.0 * x31 * y31 / l31;
    let r6 = 3.0 * x12 * y12 / l12;
    let t4 = -6.0 * x23 / l23;
    let t5 = -6.0 * x31 / l31;
    let t6 = -6.0 * x12 / l12;
    let d_l1dx = y23 / two_a;
    let d_l1dy = -x23 / two_a;
    let d_l2dx = y31 / two_a;
    let d_l2dy = -x31 / two_a;
    let pts = [[1.0 / 6.0, 1.0 / 6.0], [2.0 / 3.0, 1.0 / 6.0], [1.0 / 6.0, 2.0 / 3.0]];
    let wt = area / 3.0;
    for p in &pts {
        let l2 = p[0];
        let l1 = 1.0 - p[0] - p[1];
        // Hx, Hy interpolants (Batoz) of βx, βy ; 9 DOF [w,θx,θy]×3
        let hx = [
            p6 * (1.0 - 2.0 * l1) + (p5 - p6) * l2,
            q6 * (1.0 - 2.0 * l1) - (q5 + q6) * l2,
            -4.0 + 6.0 * (l1 + l2) + r6 * (1.0 - 2.0 * l1) - l2 * (r5 + r6),
            -p6 * (1.0 - 2.0 * l1) + (p4 + p6) * l2,
            q6 * (1.0 - 2.0 * l1) - (q6 - q4) * l2,
            -2.0 + 6.0 * l2 + r6 * (1.0 - 2.0 * l1) + l2 * (r4 - r6),
            -l2 * (p5 + p4),
            l2 * (q4 - q5),
            -l2 * (r5 - r4),
        ];
        let hy = [
            t6 * (1.0 - 2.0 * l1) + (t5 - t6) * l2,
            1.0 + r6 * (1.0 - 2.0 * l1) - l2 * (r5 + r6),
            -q6 * (1.0 - 2.0 * l1) + l2 * (q5 + q6),
            -t6 * (1.0 - 2.0 * l1) + (t4 + t6) * l2,
            -1.0 + r6 * (1.0 - 2.0 * l1) + l2 * (r4 - r6),
            -q6 * (1.0 - 2.0 * l1) + l2 * (q6 - q4),
            -l2 * (t4 + t5),
            l2 * (r4 - r5),
            l2 * (q5 - q4),
        ];
        // derivatives wrt L1, L2 (Hx, Hy linear in L)
        let dhx_dl1 = [-2.0 * p6, -2.0 * q6, 6.0 - 2.0 * r6, 2.0 * p6, -2.0 * q6, -2.0 * r6, 0.0, 0.0, 0.0];
        let dhx_dl2 = [
            p5 - p6,
            -(q5 + q6),
            6.0 - (r5 + r6),
            p4 + p6,
            -(q6 - q4),
            6.0 + (r4 - r6),
            -(p5 + p4),
            q4 - q5,
            -(r5 - r4),
        ];
        let dhy_dl1 = [-2.0 * t6, -2.0 * r6, 2.0 * q6, 2.0 * t6, -2.0 * r6, 2.0 * q6, 0.0, 0.0, 0.0];
        let dhy_dl2 = [
            t5 - t6,
            -(r5 + r6),
            q5 + q6,
            t4 + t6,
            r4 - r6,
            q6 - q4,
            -(t4 + t5),
            r4 - r5,
            q5 - q4,
        ];
        let _ = (hx, hy, d_l1dx, d_l1dy, d_l2dx, d_l2dy);
        let mut dhx_dx = [0.0; 9];
        let mut dhx_dy = [0.0; 9];
        let mut dhy_dx = [0.0; 9];
        let mut dhy_dy = [0.0; 9];
        for i in 0..9 {
            dhx_dx[i] = dhx_dl1[i] * d_l1dx + dhx_dl2[i] * d_l2dx;
            dhx_dy[i] = dhx_dl1[i] * d_l1dy + dhx_dl2[i] * d_l2dy;
            dhy_dx[i] = dhy_dl1[i] * d_l1dx + dhy_dl2[i] * d_l2dx;
            dhy_dy[i] = dhy_dl1[i] * d_l1dy + dhy_dl2[i] * d_l2dy;
        }
        // Batoz: βx ≈ Hx·d, βy ≈ Hy·d with d=[w,θx,θy]...
        // κx=βx,x  κy=βy,y  κxy=βx,y+βy,x
        // Map 9 DOF into our 18: node a → w=6a+2, θx=6a+3, θy=6a+4
        let mut b = vec![0.0; 3 * nd];
        for a in 0..3 {
            let c0 = 3 * a; // index in 9-vector
            let g = 6 * a;
            b[0 * nd + g + 2] = dhx_dx[c0];
            b[0 * nd + g + 3] = dhx_dx[c0 + 1];
            b[0 * nd + g + 4] = dhx_dx[c0 + 2];
            b[1 * nd + g + 2] = dhy_dy[c0];
            b[1 * nd + g + 3] = dhy_dy[c0 + 1];
            b[1 * nd + g + 4] = dhy_dy[c0 + 2];
            b[2 * nd + g + 2] = dhx_dy[c0] + dhy_dx[c0];
            b[2 * nd + g + 3] = dhx_dy[c0 + 1] + dhy_dx[c0 + 1];
            b[2 * nd + g + 4] = dhx_dy[c0 + 2] + dhy_dx[c0 + 2];
        }
        gemm_bt_d_b(&mut ke, nd, &b, 3, &db, wt);
    }

    // 1-pt Mindlin shear (thin-plate: small energy)
    let ds = ds_mat(e, nu, h);
    let nshp = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
    let bg = cartesian_shear_b(&nshp, &dndx, 3);
    add_shear(&mut ke, nd, &bg, &ds, area, 3);
    add_drill(&mut ke, nd, 3, e, h, area);
    Ok((ke, area))
}

fn s6_local(xy: &[[f64; 2]], e: f64, nu: f64, h: f64) -> Result<(Vec<f64>, f64)> {
    let nn = 6usize;
    let nd = 36usize;
    let mut ke = vec![0.0; nd * nd];
    let dm0 = d_plane_stress(e, nu)?;
    let mut dm = [0.0; 9];
    let mut db = [0.0; 9];
    for i in 0..9 {
        dm[i] = dm0[i] * h;
        db[i] = dm0[i] * h * h * h / 12.0;
    }
    let ds = ds_mat(e, nu, h);
    let mut area = 0.0;
    let pts = [[1.0 / 6.0, 1.0 / 6.0], [2.0 / 3.0, 1.0 / 6.0], [1.0 / 6.0, 2.0 / 3.0]];
    let w = 1.0 / 6.0;
    for p in &pts {
        let (nshp, dn) = tri6_shape(p[0], p[1]);
        let mut dna = [[0.0; 2]; 6];
        dna.copy_from_slice(&dn);
        let (_, det, dndx) = jac_xy(xy, &dna, nn)?;
        if det <= 0.0 {
            return err("S6: negative Jakobideterminante.");
        }
        add_membrane(&mut ke, nd, nn, &dndx, &dm, w * det);
        add_bending(&mut ke, nd, nn, &dndx, &db, w * det);
        let bg = cartesian_shear_b(&nshp, &dndx, nn);
        add_shear(&mut ke, nd, &bg, &ds, w * det, nn);
        area += w * det;
    }
    add_drill(&mut ke, nd, nn, e, h, area);
    Ok((ke, area))
}

pub fn stiffness(kind: ElemKind, xyz: &[[f64; 3]], e: f64, nu: f64, h: f64) -> Result<(Vec<f64>, f64)> {
    let nn = kind.nnodes();
    if xyz.len() < nn {
        return err("Schale: zu wenige Knoten.");
    }
    if h <= 0.0 {
        return err("SHELL SECTION: Dicke muss positiv sein.");
    }
    let (e1, e2, e3) = local_frame(xyz, nn)?;
    let xy = project_xy(xyz, e1, e2, nn);
    let (mut ke, area) = match kind {
        ElemKind::Shell4 | ElemKind::Shell4R => {
            let mut p = [[0.0; 2]; 4];
            p.copy_from_slice(&xy[..4]);
            s4_local(&p, e, nu, h)?
        }
        ElemKind::Shell8 | ElemKind::Shell8R => s8_local(&xy, e, nu, h, kind.reduced_int())?,
        ElemKind::Shell3 => {
            let mut p = [[0.0; 2]; 3];
            p.copy_from_slice(&xy[..3]);
            s3_local(&p, e, nu, h)?
        }
        ElemKind::Shell6 => s6_local(&xy, e, nu, h)?,
        _ => return err("Kein Schalenelement."),
    };
    rotate_ke(&mut ke, nn, e1, e2, e3);
    let _ = e3;
    Ok((ke, area))
}

fn tensor_rotate(sl: [f64; 6], e1: [f64; 3], e2: [f64; 3], e3: [f64; 3]) -> [f64; 6] {
    // local Voigt [sxx, syy, szz, sxy, syz, szx] with szz=0 typically
    let mut sm = [[0.0; 3]; 3];
    sm[0][0] = sl[0];
    sm[1][1] = sl[1];
    sm[2][2] = sl[2];
    sm[0][1] = sl[3];
    sm[1][0] = sl[3];
    sm[1][2] = sl[4];
    sm[2][1] = sl[4];
    sm[0][2] = sl[5];
    sm[2][0] = sl[5];
    let r = [e1, e2, e3]; // rows are local axes in global
    // σg = R^T σl R  where R columns = e1,e2,e3 so R[i][j] = e_j[i]
    let mut sg = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    sg[i][j] += r[k][i] * sm[k][l] * r[l][j];
                }
            }
        }
    }
    [sg[0][0], sg[1][1], sg[2][2], sg[0][1], sg[1][2], sg[0][2]]
}

fn local_fiber_stress(
    dndx: &[[f64; 2]],
    nshp: &[f64],
    ul: &[f64],
    dm0: &[f64],
    h: f64,
    nn: usize,
    e: f64,
    nu: f64,
) -> [f64; 6] {
    // membrane strain
    let mut eps = [0.0; 3];
    let mut kap = [0.0; 3];
    let mut gam = [0.0; 2];
    for a in 0..nn {
        let u = ul[6 * a];
        let v = ul[6 * a + 1];
        let w = ul[6 * a + 2];
        let tx = ul[6 * a + 3];
        let ty = ul[6 * a + 4];
        eps[0] += dndx[a][0] * u;
        eps[1] += dndx[a][1] * v;
        eps[2] += dndx[a][1] * u + dndx[a][0] * v;
        kap[0] += dndx[a][0] * ty;
        kap[1] += -dndx[a][1] * tx;
        kap[2] += dndx[a][1] * ty - dndx[a][0] * tx;
        gam[0] += dndx[a][0] * w + nshp[a] * ty;
        gam[1] += dndx[a][1] * w - nshp[a] * tx;
    }
    let z = 0.5 * h;
    let mut em = [0.0; 3];
    for i in 0..3 {
        em[i] = eps[i] + z * kap[i];
    }
    let mut sm = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            sm[i] += dm0[i * 3 + j] * em[j];
        }
    }
    let g = e / (2.0 * (1.0 + nu)) * K_SHEAR;
    [sm[0], sm[1], 0.0, sm[2], g * gam[1], g * gam[0]]
}

pub fn nodal_stress(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    h: f64,
) -> Result<Vec<[f64; 6]>> {
    let nn = kind.nnodes();
    let (e1, e2, e3) = local_frame(xyz, nn)?;
    let xy = project_xy(xyz, e1, e2, nn);
    let ul = rotate_vec_to_local(ue, nn, e1, e2, e3);
    let dm0 = d_plane_stress(e, nu)?;
    let mut out = vec![[0.0; 6]; nn];
    match kind {
        ElemKind::Shell4 | ElemKind::Shell4R => {
            for a in 0..4 {
                let (nshp, dn) = quad4_shape(QUAD_XI[a][0], QUAD_XI[a][1]);
                let mut p = [[0.0; 2]; 4];
                p.copy_from_slice(&xy[..4]);
                let (_, _, dndx) = jac_xy(&p, &dn, 4)?;
                let sl = local_fiber_stress(&dndx, &nshp, &ul, &dm0, h, 4, e, nu);
                out[a] = tensor_rotate(sl, e1, e2, e3);
            }
        }
        ElemKind::Shell8 | ElemKind::Shell8R => {
            const Q8: [[f64; 2]; 8] = [
                [-1.0, -1.0],
                [1.0, -1.0],
                [1.0, 1.0],
                [-1.0, 1.0],
                [0.0, -1.0],
                [1.0, 0.0],
                [0.0, 1.0],
                [-1.0, 0.0],
            ];
            for a in 0..8 {
                let (nshp, dn) = quad8_shape(Q8[a][0], Q8[a][1]);
                let mut dna = [[0.0; 2]; 8];
                dna.copy_from_slice(&dn);
                let (_, _, dndx) = jac_xy(&xy, &dna, 8)?;
                let sl = local_fiber_stress(&dndx, &nshp, &ul, &dm0, h, 8, e, nu);
                out[a] = tensor_rotate(sl, e1, e2, e3);
            }
        }
        ElemKind::Shell3 => {
            let mut p = [[0.0; 2]; 3];
            p.copy_from_slice(&xy[..3]);
            let x1 = p[0][0];
            let y1 = p[0][1];
            let x2 = p[1][0];
            let y2 = p[1][1];
            let x3 = p[2][0];
            let y3 = p[2][1];
            let two_a = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
            let mut dndx = [[0.0; 2]; 3];
            dndx[0] = [(y2 - y3) / two_a, (x3 - x2) / two_a];
            dndx[1] = [(y3 - y1) / two_a, (x1 - x3) / two_a];
            dndx[2] = [(y1 - y2) / two_a, (x2 - x1) / two_a];
            let nshp = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
            let sl = local_fiber_stress(&dndx, &nshp, &ul, &dm0, h, 3, e, nu);
            for a in 0..3 {
                out[a] = tensor_rotate(sl, e1, e2, e3);
            }
        }
        ElemKind::Shell6 => {
            let rst = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.5, 0.0], [0.5, 0.5], [0.0, 0.5]];
            for a in 0..6 {
                let (nshp, dn) = tri6_shape(rst[a][0], rst[a][1]);
                let mut dna = [[0.0; 2]; 6];
                dna.copy_from_slice(&dn);
                let (_, _, dndx) = jac_xy(&xy, &dna, 6)?;
                let sl = local_fiber_stress(&dndx, &nshp, &ul, &dm0, h, 6, e, nu);
                out[a] = tensor_rotate(sl, e1, e2, e3);
            }
        }
        _ => return err("Kein Schalenelement."),
    }
    Ok(out)
}

pub fn pressure_force(kind: ElemKind, xyz: &[[f64; 3]], p: f64) -> Result<Vec<f64>> {
    let nn = kind.nnodes();
    let (e1, e2, e3) = local_frame(xyz, nn)?;
    let xy = project_xy(xyz, e1, e2, nn);
    let mut fe = vec![0.0; 6 * nn];
    let apply = |nshp: &[f64], det: f64, w: f64, fe: &mut [f64], e3: [f64; 3], nn: usize| {
        let f = p * w * det;
        for a in 0..nn {
            fe[6 * a] += nshp[a] * f * e3[0];
            fe[6 * a + 1] += nshp[a] * f * e3[1];
            fe[6 * a + 2] += nshp[a] * f * e3[2];
        }
    };
    match kind {
        ElemKind::Shell4 | ElemKind::Shell4R => {
            let mut pxy = [[0.0; 2]; 4];
            pxy.copy_from_slice(&xy[..4]);
            for &xi in &[-G2, G2] {
                for &eta in &[-G2, G2] {
                    let (nshp, dn) = quad4_shape(xi, eta);
                    let (_, det, _) = jac_xy(&pxy, &dn, 4)?;
                    apply(&nshp, det, 1.0, &mut fe, e3, 4);
                }
            }
        }
        ElemKind::Shell8 | ElemKind::Shell8R => {
            for i in 0..3 {
                for j in 0..3 {
                    let (nshp, dn) = quad8_shape(G3[i], G3[j]);
                    let mut dna = [[0.0; 2]; 8];
                    dna.copy_from_slice(&dn);
                    let (_, det, _) = jac_xy(&xy, &dna, 8)?;
                    apply(&nshp, det, W3[i] * W3[j], &mut fe, e3, 8);
                }
            }
        }
        ElemKind::Shell3 => {
            let (nshp, _) = tri3_shape();
            let mut pxy = [[0.0; 2]; 3];
            pxy.copy_from_slice(&xy[..3]);
            let two_a = pxy[0][0] * (pxy[1][1] - pxy[2][1])
                + pxy[1][0] * (pxy[2][1] - pxy[0][1])
                + pxy[2][0] * (pxy[0][1] - pxy[1][1]);
            apply(&nshp, two_a * 0.5, 1.0, &mut fe, e3, 3);
        }
        ElemKind::Shell6 => {
            let pts = [[1.0 / 6.0, 1.0 / 6.0], [2.0 / 3.0, 1.0 / 6.0], [1.0 / 6.0, 2.0 / 3.0]];
            for gp in &pts {
                let (nshp, dn) = tri6_shape(gp[0], gp[1]);
                let mut dna = [[0.0; 2]; 6];
                dna.copy_from_slice(&dn);
                let (_, det, _) = jac_xy(&xy, &dna, 6)?;
                apply(&nshp, det, 1.0 / 6.0, &mut fe, e3, 6);
            }
        }
        _ => {}
    }
    let _ = (e1, e2);
    Ok(fe)
}

pub fn body_force(kind: ElemKind, xyz: &[[f64; 3]], bx: f64, by: f64, bz: f64, h: f64) -> Result<Vec<f64>> {
    let nn = kind.nnodes();
    let (e1, e2, e3) = local_frame(xyz, nn)?;
    let xy = project_xy(xyz, e1, e2, nn);
    let mut fe = vec![0.0; 6 * nn];
    let g = [bx, by, bz];
    let apply = |nshp: &[f64], det: f64, w: f64, fe: &mut [f64], nn: usize| {
        let f = h * w * det;
        for a in 0..nn {
            fe[6 * a] += nshp[a] * f * g[0];
            fe[6 * a + 1] += nshp[a] * f * g[1];
            fe[6 * a + 2] += nshp[a] * f * g[2];
        }
    };
    match kind {
        ElemKind::Shell4 | ElemKind::Shell4R => {
            let mut pxy = [[0.0; 2]; 4];
            pxy.copy_from_slice(&xy[..4]);
            for &xi in &[-G2, G2] {
                for &eta in &[-G2, G2] {
                    let (nshp, dn) = quad4_shape(xi, eta);
                    let (_, det, _) = jac_xy(&pxy, &dn, 4)?;
                    apply(&nshp, det, 1.0, &mut fe, 4);
                }
            }
        }
        ElemKind::Shell8 | ElemKind::Shell8R => {
            for i in 0..3 {
                for j in 0..3 {
                    let (nshp, dn) = quad8_shape(G3[i], G3[j]);
                    let mut dna = [[0.0; 2]; 8];
                    dna.copy_from_slice(&dn);
                    let (_, det, _) = jac_xy(&xy, &dna, 8)?;
                    apply(&nshp, det, W3[i] * W3[j], &mut fe, 8);
                }
            }
        }
        ElemKind::Shell3 => {
            let (nshp, _) = tri3_shape();
            let mut pxy = [[0.0; 2]; 3];
            pxy.copy_from_slice(&xy[..3]);
            let two_a = pxy[0][0] * (pxy[1][1] - pxy[2][1])
                + pxy[1][0] * (pxy[2][1] - pxy[0][1])
                + pxy[2][0] * (pxy[0][1] - pxy[1][1]);
            apply(&nshp, two_a.abs() * 0.5, 1.0, &mut fe, 3);
        }
        ElemKind::Shell6 => {
            let pts = [[1.0 / 6.0, 1.0 / 6.0], [2.0 / 3.0, 1.0 / 6.0], [1.0 / 6.0, 2.0 / 3.0]];
            for gp in &pts {
                let (nshp, dn) = tri6_shape(gp[0], gp[1]);
                let mut dna = [[0.0; 2]; 6];
                dna.copy_from_slice(&dn);
                let (_, det, _) = jac_xy(&xy, &dna, 6)?;
                apply(&nshp, det, 1.0 / 6.0, &mut fe, 6);
            }
        }
        _ => {}
    }
    let _ = e3;
    Ok(fe)
}

/// Membrane elements M3D3/M3D4/M3D4R/M3D6/M3D8: plane-stress in the local tangent,
/// 3 translational DOF, no bending.
pub fn membrane_stiffness(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    e: f64,
    nu: f64,
    h: f64,
) -> Result<(Vec<f64>, f64)> {
    let nn = kind.nnodes();
    if xyz.len() < nn {
        return err("Membran: zu wenige Knoten.");
    }
    if h <= 0.0 {
        return err("MEMBRANE SECTION: Dicke muss positiv sein.");
    }
    let (e1, e2, _e3) = local_frame(xyz, nn)?;
    let xy = project_xy(xyz, e1, e2, nn);
    let (ke2, area) = match kind {
        ElemKind::Mem4 | ElemKind::Mem4R => {
            let mut p = [[0.0; 2]; 4];
            p.copy_from_slice(&xy[..4]);
            mem_quad4(&p, e, nu, h, kind.reduced_int())?
        }
        ElemKind::Mem8 => {
            let mut p = [[0.0; 2]; 8];
            p.copy_from_slice(&xy[..8]);
            crate::quadratic::quad8_stiffness(&p, e, nu, h, false, false)?
        }
        ElemKind::Mem3 => {
            let mut p = [[0.0; 2]; 3];
            p.copy_from_slice(&xy[..3]);
            mem_tri3(&p, e, nu, h)?
        }
        ElemKind::Mem6 => {
            let mut p = [[0.0; 2]; 6];
            p.copy_from_slice(&xy[..6]);
            crate::quadratic::tri6_stiffness(&p, e, nu, h, false)?
        }
        _ => return err("Kein Membranelement."),
    };
    let nd = 3 * nn;
    let mut ke = vec![0.0; nd * nd];
    let n2 = 2 * nn;
    for a in 0..nn {
        for b in 0..nn {
            for i in 0..2 {
                for j in 0..2 {
                    let v = ke2[(2 * a + i) * n2 + (2 * b + j)];
                    // T columns are e1, e2
                    let ti = if i == 0 { e1 } else { e2 };
                    let tj = if j == 0 { e1 } else { e2 };
                    for p in 0..3 {
                        for q in 0..3 {
                            ke[(3 * a + p) * nd + (3 * b + q)] += ti[p] * v * tj[q];
                        }
                    }
                }
            }
        }
    }
    Ok((ke, area))
}

fn mem_quad4(xy: &[[f64; 2]; 4], e: f64, nu: f64, h: f64, reduced: bool) -> Result<(Vec<f64>, f64)> {
    let d = d_plane_stress(e, nu)?;
    let n = 8usize;
    let mut ke = vec![0.0; n * n];
    let mut area = 0.0;
    let gps: Vec<(f64, f64, f64)> = if reduced {
        vec![(0.0, 0.0, 4.0)]
    } else {
        let mut o = Vec::new();
        for &xi in &[-G2, G2] {
            for &eta in &[-G2, G2] {
                o.push((xi, eta, 1.0));
            }
        }
        o
    };
    for (xi, eta, w0) in gps {
        let (nshp, dn) = quad4_shape(xi, eta);
        let _ = nshp;
        let (_, det, dndx) = jac_xy(&xy.to_vec(), &dn, 4)?;
        if det <= 0.0 {
            return err("M3D4: negative Jakobideterminante.");
        }
        let mut b = vec![0.0; 3 * n];
        let mut d2 = [[0.0; 2]; 4];
        for i in 0..4 {
            d2[i] = dndx[i];
        }
        fill_b2(&mut b, 4, &d2);
        gemm_bt_d_b(&mut ke, n, &b, 3, &d, h * w0 * det);
        area += w0 * det;
    }
    Ok((ke, area))
}

fn mem_tri3(xy: &[[f64; 2]; 3], e: f64, nu: f64, h: f64) -> Result<(Vec<f64>, f64)> {
    let x1 = xy[0][0];
    let y1 = xy[0][1];
    let x2 = xy[1][0];
    let y2 = xy[1][1];
    let x3 = xy[2][0];
    let y3 = xy[2][1];
    let two_a = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
    if two_a <= 0.0 {
        return err("M3D3: nicht-positive Fläche.");
    }
    let a = 0.5 * two_a;
    let mut dndx = [[0.0; 2]; 3];
    dndx[0] = [(y2 - y3) / two_a, (x3 - x2) / two_a];
    dndx[1] = [(y3 - y1) / two_a, (x1 - x3) / two_a];
    dndx[2] = [(y1 - y2) / two_a, (x2 - x1) / two_a];
    let d = d_plane_stress(e, nu)?;
    let n = 6usize;
    let mut ke = vec![0.0; n * n];
    let mut b = vec![0.0; 3 * n];
    fill_b2(&mut b, 3, &dndx);
    gemm_bt_d_b(&mut ke, n, &b, 3, &d, h * a);
    Ok((ke, a))
}

pub fn membrane_nodal_stress(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    h: f64,
) -> Result<Vec<[f64; 6]>> {
    let _ = h;
    let nn = kind.nnodes();
    let (e1, e2, _e3) = local_frame(xyz, nn)?;
    let xy = project_xy(xyz, e1, e2, nn);
    let mut u2 = vec![0.0; 2 * nn];
    for a in 0..nn {
        u2[2 * a] = ue[3 * a] * e1[0] + ue[3 * a + 1] * e1[1] + ue[3 * a + 2] * e1[2];
        u2[2 * a + 1] = ue[3 * a] * e2[0] + ue[3 * a + 1] * e2[1] + ue[3 * a + 2] * e2[2];
    }
    let d = d_plane_stress(e, nu)?;
    let n2 = 2 * nn;
    let mut slocal = vec![[0.0; 3]; nn];
    for a in 0..nn {
        let (dndx, ok) = match kind {
            ElemKind::Mem4 | ElemKind::Mem4R => {
                let (nshp, dn) = quad4_shape(QUAD_XI[a.min(3)][0], QUAD_XI[a.min(3)][1]);
                let _ = nshp;
                let r = jac_xy(&xy, &dn, 4);
                match r {
                    Ok((_, det, dndx)) if det > 0.0 => (dndx, true),
                    _ => (vec![[0.0; 2]; nn], false),
                }
            }
            ElemKind::Mem8 => {
                let (nshp, dn) = quad8_shape(
                    crate::quadratic::QUAD8_XI[a.min(7)][0],
                    crate::quadratic::QUAD8_XI[a.min(7)][1],
                );
                let _ = nshp;
                match jac_xy(&xy, &dn, 8) {
                    Ok((_, det, dndx)) if det > 0.0 => (dndx, true),
                    _ => (vec![[0.0; 2]; nn], false),
                }
            }
            ElemKind::Mem3 => {
                let two_a = xy[0][0] * (xy[1][1] - xy[2][1])
                    + xy[1][0] * (xy[2][1] - xy[0][1])
                    + xy[2][0] * (xy[0][1] - xy[1][1]);
                if two_a <= 0.0 {
                    (vec![[0.0; 2]; 3], false)
                } else {
                    (
                        vec![
                            [(xy[1][1] - xy[2][1]) / two_a, (xy[2][0] - xy[1][0]) / two_a],
                            [(xy[2][1] - xy[0][1]) / two_a, (xy[0][0] - xy[2][0]) / two_a],
                            [(xy[0][1] - xy[1][1]) / two_a, (xy[1][0] - xy[0][0]) / two_a],
                        ],
                        true,
                    )
                }
            }
            ElemKind::Mem6 => {
                let xi = if a < 3 {
                    [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]][a]
                } else {
                    [[0.5, 0.5], [0.0, 0.5], [0.5, 0.0]][a - 3]
                };
                let (nshp, dn) = crate::quadratic::tri6_shape(xi[0], xi[1]);
                let _ = nshp;
                match jac_xy(&xy, &dn, 6) {
                    Ok((_, det, dndx)) if det > 0.0 => (dndx, true),
                    _ => (vec![[0.0; 2]; nn], false),
                }
            }
            _ => (vec![[0.0; 2]; nn], false),
        };
        if !ok {
            continue;
        }
        let mut b = vec![0.0; 3 * n2];
        let mut d2 = vec![[0.0; 2]; nn];
        for i in 0..nn {
            d2[i] = dndx[i];
        }
        fill_b2(&mut b, nn, &d2);
        let s = crate::elem::sigma_from_b(&b, 3, n2, &d, &u2);
        slocal[a] = [s[0], s[1], s[2]];
    }
    let mut out = vec![[0.0; 6]; nn];
    for a in 0..nn {
        let s11 = slocal[a][0];
        let s22 = slocal[a][1];
        let s12 = slocal[a][2];
        let t = |i: usize, j: usize| {
            s11 * e1[i] * e1[j] + s22 * e2[i] * e2[j] + s12 * (e1[i] * e2[j] + e2[i] * e1[j])
        };
        out[a] = [t(0, 0), t(1, 1), t(2, 2), t(0, 1), t(1, 2), t(2, 0)];
    }
    Ok(out)
}

pub fn membrane_pressure(kind: ElemKind, xyz: &[[f64; 3]], p: f64) -> Result<Vec<f64>> {
    let nn = kind.nnodes();
    let mut fe6 = pressure_force(
        match kind {
            ElemKind::Mem4 | ElemKind::Mem4R => ElemKind::Shell4,
            ElemKind::Mem8 => ElemKind::Shell8,
            ElemKind::Mem3 => ElemKind::Shell3,
            ElemKind::Mem6 => ElemKind::Shell6,
            k => k,
        },
        xyz,
        p,
    )?;
    let mut fe = vec![0.0; 3 * nn];
    for a in 0..nn {
        fe[3 * a] = fe6[6 * a];
        fe[3 * a + 1] = fe6[6 * a + 1];
        fe[3 * a + 2] = fe6[6 * a + 2];
    }
    Ok(fe)
}

/// Co-rotational shell: linear stiffness on the current mid-surface,
/// internal force from deformational displacement (rigid motion removed).
pub fn stiffness_nl(
    kind: ElemKind,
    xyz0: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    h: f64,
) -> Result<(Vec<f64>, Vec<f64>, [f64; 6])> {
    let nn = kind.nnodes();
    let nd = 6 * nn;
    let mut xyz = vec![[0.0; 3]; nn];
    for a in 0..nn.min(xyz0.len()) {
        xyz[a] = [
            xyz0[a][0] + ue.get(6 * a).copied().unwrap_or(0.0),
            xyz0[a][1] + ue.get(6 * a + 1).copied().unwrap_or(0.0),
            xyz0[a][2] + ue.get(6 * a + 2).copied().unwrap_or(0.0),
        ];
    }
    let (e1, e2, e3) = local_frame(&xyz, nn)?;
    let (e10, e20, e30) = local_frame(xyz0, nn)?;
    let mut c0 = [0.0; 3];
    let mut c1 = [0.0; 3];
    for a in 0..nn {
        for d in 0..3 {
            c0[d] += xyz0[a][d];
            c1[d] += xyz[a][d];
        }
    }
    let invn = 1.0 / nn as f64;
    for d in 0..3 {
        c0[d] *= invn;
        c1[d] *= invn;
    }
    let mut r = [[0.0; 3]; 3];
    let t0 = [e10, e20, e30];
    let t1 = [e1, e2, e3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = t1[0][i] * t0[0][j] + t1[1][i] * t0[1][j] + t1[2][i] * t0[2][j];
        }
    }
    let th = rotvec_from_r_shell(r);
    let mut udef = vec![0.0; nd];
    for a in 0..nn {
        let d0 = [
            xyz0[a][0] - c0[0],
            xyz0[a][1] - c0[1],
            xyz0[a][2] - c0[2],
        ];
        let rd = [
            r[0][0] * d0[0] + r[0][1] * d0[1] + r[0][2] * d0[2],
            r[1][0] * d0[0] + r[1][1] * d0[1] + r[1][2] * d0[2],
            r[2][0] * d0[0] + r[2][1] * d0[1] + r[2][2] * d0[2],
        ];
        let urig = [
            c1[0] + rd[0] - xyz0[a][0],
            c1[1] + rd[1] - xyz0[a][1],
            c1[2] + rd[2] - xyz0[a][2],
        ];
        udef[6 * a] = ue.get(6 * a).copied().unwrap_or(0.0) - urig[0];
        udef[6 * a + 1] = ue.get(6 * a + 1).copied().unwrap_or(0.0) - urig[1];
        udef[6 * a + 2] = ue.get(6 * a + 2).copied().unwrap_or(0.0) - urig[2];
        udef[6 * a + 3] = ue.get(6 * a + 3).copied().unwrap_or(0.0) - th[0];
        udef[6 * a + 4] = ue.get(6 * a + 4).copied().unwrap_or(0.0) - th[1];
        udef[6 * a + 5] = ue.get(6 * a + 5).copied().unwrap_or(0.0) - th[2];
    }
    let (ke0, _) = stiffness(kind, xyz0, e, nu, h)?;
    let (ke, _) = stiffness(kind, &xyz, e, nu, h)?;
    let mut fe = vec![0.0; nd];
    for i in 0..nd {
        let mut s = 0.0;
        for j in 0..nd {
            s += ke0.get(i * nd + j).copied().unwrap_or(0.0) * udef[j];
        }
        fe[i] = s;
    }
    let stress = nodal_stress(kind, &xyz, &udef, e, nu, h)
        .ok()
        .and_then(|v| v.first().copied())
        .unwrap_or([0.0; 6]);
    let _ = (e2, e20);
    Ok((ke, fe, stress))
}

fn rotvec_from_r_shell(r: [[f64; 3]; 3]) -> [f64; 3] {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn s4_pressure_integrates_to_area() {
        let xyz = [
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.0, 10.0, 0.0],
            [0.0, 10.0, 0.0],
        ];
        let fe = pressure_force(ElemKind::Shell4, &xyz, 1.0).unwrap();
        let mut fz = 0.0;
        for a in 0..4 {
            fz += fe[6 * a + 2];
        }
        assert!((fz - 100.0).abs() < 1e-6, "Fz={fz}");
    }
}
