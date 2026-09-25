//! Curved MITC4 shell (Dvorkin–Bathe).
//!
//! The midsurface is the bilinear patch. Directors are the nodal surface normals
//! (averaged by the assembler). Membrane and curvature strains use the surface
//! metric, so a faceted sphere still carries the inextensional bending mode.
//! Transverse shear is tied at the four edge midpoints.

use crate::elem::{d_plane_stress, gemm_bt_d_b, invert2, G2};
use crate::error::{err, Result};
use crate::ortho::ShellLaw;

const K_SHEAR: f64 = 5.0 / 6.0;
const DRILL: f64 = 1.0e-8;
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
        return err("S4: Normale ist degeneriert.");
    }
    Ok(scale(a, 1.0 / n))
}
fn unit3(d: usize) -> [f64; 3] {
    match d {
        0 => [1.0, 0.0, 0.0],
        1 => [0.0, 1.0, 0.0],
        _ => [0.0, 0.0, 1.0],
    }
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

fn corner_normals(xyz: &[[f64; 3]; 4]) -> Result<[[f64; 3]; 4]> {
    let mut n = [[0.0; 3]; 4];
    for a in 0..4 {
        let (_, dn) = quad4_shape(QUAD_XI[a][0], QUAD_XI[a][1]);
        let mut gxi = [0.0; 3];
        let mut geta = [0.0; 3];
        for b in 0..4 {
            gxi = add(gxi, scale(xyz[b], dn[b][0]));
            geta = add(geta, scale(xyz[b], dn[b][1]));
        }
        n[a] = normalize(cross(gxi, geta))?;
    }
    Ok(n)
}

fn resolve_directors(xyz: &[[f64; 3]; 4], directors: Option<&[[f64; 3]]>) -> Result<[[f64; 3]; 4]> {
    let fb = corner_normals(xyz)?;
    let Some(d) = directors else {
        return Ok(fb);
    };
    if d.len() < 4 {
        return Ok(fb);
    }
    let mut an = [[0.0; 3]; 4];
    for a in 0..4 {
        an[a] = if norm(d[a]) > 1e-8 {
            normalize(d[a])?
        } else {
            fb[a]
        };
    }
    Ok(an)
}

struct Kin {
    gxi: [f64; 3],
    geta: [f64; 3],
    v: [f64; 3],
    vxi: [f64; 3],
    veta: [f64; 3],
    nshp: [f64; 4],
    dn: [[f64; 2]; 4],
    lambda: f64,
    e1: [f64; 3],
    e2: [f64; 3],
    e3: [f64; 3],
    xix: f64,
    xiy: f64,
    etx: f64,
    ety: f64,
    da: f64,
}

fn kin(xyz: &[[f64; 3]; 4], an: &[[f64; 3]; 4], xi: f64, eta: f64) -> Result<Kin> {
    let (nshp, dn) = quad4_shape(xi, eta);
    let mut gxi = [0.0; 3];
    let mut geta = [0.0; 3];
    let mut v = [0.0; 3];
    let mut vxi = [0.0; 3];
    let mut veta = [0.0; 3];
    for a in 0..4 {
        gxi = add(gxi, scale(xyz[a], dn[a][0]));
        geta = add(geta, scale(xyz[a], dn[a][1]));
        v = add(v, scale(an[a], nshp[a]));
        vxi = add(vxi, scale(an[a], dn[a][0]));
        veta = add(veta, scale(an[a], dn[a][1]));
    }
    let lambda = norm(v);
    if lambda < 1e-14 {
        return err("S4: Direktor verschwindet.");
    }
    let e3 = scale(v, 1.0 / lambda);
    let t1 = sub(gxi, scale(e3, dot(gxi, e3)));
    if norm(t1) < 1e-14 {
        return err("S4: entartete Tangente.");
    }
    let e1 = normalize(t1)?;
    let e2 = cross(e3, e1);
    let j = [
        [dot(gxi, e1), dot(geta, e1)],
        [dot(gxi, e2), dot(geta, e2)],
    ];
    let (inv, det) = invert2(j)?;
    if det <= 0.0 {
        return err("S4: negative Jakobideterminante.");
    }
    Ok(Kin {
        gxi,
        geta,
        v,
        vxi,
        veta,
        nshp,
        dn,
        lambda,
        e1,
        e2,
        e3,
        xix: inv[0][0],
        xiy: inv[0][1],
        etx: inv[1][0],
        ety: inv[1][1],
        da: norm(cross(gxi, geta)),
    })
}

fn cov_to_cart(c: &[f64], k: &Kin) -> Vec<f64> {
    let nd = 24usize;
    let mut b = vec![0.0; 3 * nd];
    let (xx, xy, yx, yy) = (k.xix, k.xiy, k.etx, k.ety);
    for col in 0..nd {
        let exx = c[col];
        let eyy = c[nd + col];
        let exy = c[2 * nd + col];
        b[col] = xx * xx * exx + yx * yx * eyy + xx * yx * exy;
        b[nd + col] = xy * xy * exx + yy * yy * eyy + xy * yy * exy;
        b[2 * nd + col] = 2.0 * xx * xy * exx + 2.0 * yx * yy * eyy + (xx * yy + yx * xy) * exy;
    }
    b
}

fn mb(k: &Kin, an: &[[f64; 3]; 4]) -> (Vec<f64>, Vec<f64>) {
    let nd = 24usize;
    let mut cm = vec![0.0; 3 * nd];
    let mut ck = vec![0.0; 3 * nd];
    for a in 0..4 {
        let nx = k.dn[a][0];
        let ny = k.dn[a][1];
        for d in 0..3 {
            let ed = unit3(d);
            let col = 6 * a + d;
            let gxi_d = dot(k.gxi, ed);
            let get_d = dot(k.geta, ed);
            cm[col] += nx * gxi_d;
            cm[nd + col] += ny * get_d;
            cm[2 * nd + col] += ny * gxi_d + nx * get_d;
            ck[col] += nx * dot(k.vxi, ed);
            ck[nd + col] += ny * dot(k.veta, ed);
            ck[2 * nd + col] += ny * dot(k.vxi, ed) + nx * dot(k.veta, ed);
        }
        let n_x_gxi = cross(an[a], k.gxi);
        let n_x_get = cross(an[a], k.geta);
        for d in 0..3 {
            let ed = unit3(d);
            let col = 6 * a + 3 + d;
            let cxi = dot(n_x_gxi, ed);
            let cet = dot(n_x_get, ed);
            ck[col] += nx * cxi;
            ck[nd + col] += ny * cet;
            ck[2 * nd + col] += ny * cxi + nx * cet;
        }
    }
    let lam = k.lambda;
    for v in &mut ck {
        *v /= lam;
    }
    (cov_to_cart(&cm, k), cov_to_cart(&ck, k))
}

fn shear_row(k: &Kin, an: &[[f64; 3]; 4], eta_comp: bool) -> [f64; 24] {
    let mut row = [0.0; 24];
    let (g, comp) = if eta_comp {
        (k.geta, 1usize)
    } else {
        (k.gxi, 0usize)
    };
    let lam = k.lambda;
    for a in 0..4 {
        let nxg = cross(an[a], g);
        for d in 0..3 {
            let ed = unit3(d);
            row[6 * a + d] += k.dn[a][comp] * dot(k.v, ed) / lam;
            row[6 * a + 3 + d] += k.nshp[a] * dot(nxg, ed) / lam;
        }
    }
    row
}

fn shear_b(k: &Kin, xi: f64, eta: f64, tied: &[[f64; 24]; 4]) -> Vec<f64> {
    let mut gs = [0.0; 24];
    let mut gt = [0.0; 24];
    for i in 0..24 {
        gs[i] = 0.5 * (1.0 + eta) * tied[0][i] + 0.5 * (1.0 - eta) * tied[1][i];
        gt[i] = 0.5 * (1.0 + xi) * tied[2][i] + 0.5 * (1.0 - xi) * tied[3][i];
    }
    let mut b = vec![0.0; 48];
    for col in 0..24 {
        b[col] = k.xix * gs[col] + k.etx * gt[col];
        b[24 + col] = k.xiy * gs[col] + k.ety * gt[col];
    }
    b
}

fn tied_rows(xyz: &[[f64; 3]; 4], an: &[[f64; 3]; 4]) -> Result<[[f64; 24]; 4]> {
    Ok([
        shear_row(&kin(xyz, an, 0.0, 1.0)?, an, false),
        shear_row(&kin(xyz, an, 0.0, -1.0)?, an, false),
        shear_row(&kin(xyz, an, 1.0, 0.0)?, an, true),
        shear_row(&kin(xyz, an, -1.0, 0.0)?, an, true),
    ])
}

fn tensor_rotate(sl: [f64; 6], e1: [f64; 3], e2: [f64; 3], e3: [f64; 3]) -> [f64; 6] {
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
    let r = [e1, e2, e3];
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

pub fn s4_ke(
    xyz: &[[f64; 3]],
    e: f64,
    nu: f64,
    h: f64,
    directors: Option<&[[f64; 3]]>,
) -> Result<(Vec<f64>, f64)> {
    s4_ke_bi(xyz, e, nu, h, directors, 1.0)
}

/// `bi` is Nastran PSHELL 12I/T³. 1 recovers the homogeneous plate.
pub fn s4_ke_bi(
    xyz: &[[f64; 3]],
    e: f64,
    nu: f64,
    h: f64,
    directors: Option<&[[f64; 3]]>,
    bi: f64,
) -> Result<(Vec<f64>, f64)> {
    let law = crate::ortho::isotropic_shell(e, nu, h, bi)?;
    s4_ke_law(xyz, &law, directors, 1.0)
}

/// MITC4 with a pre-integrated plate law (MAT1, MAT2, MAT8 or PCOMP).
/// `k6rot` multiplies the drilling stabilizer (PARAM K6ROT, default 1).
pub fn s4_ke_law(
    xyz: &[[f64; 3]],
    law: &ShellLaw,
    directors: Option<&[[f64; 3]]>,
    k6rot: f64,
) -> Result<(Vec<f64>, f64)> {
    if xyz.len() < 4 {
        return err("S4: zu wenige Knoten.");
    }
    let mut p = [[0.0; 3]; 4];
    p.copy_from_slice(&xyz[..4]);
    let an = resolve_directors(&p, directors)?;
    let nd = 24usize;
    let mut ke = vec![0.0; nd * nd];
    let tied = tied_rows(&p, &an)?;
    let mut area = 0.0;
    let mut e3_avg = [0.0; 3];
    for &xi in &[-G2, G2] {
        for &eta in &[-G2, G2] {
            let k = kin(&p, &an, xi, eta)?;
            let (bm, bb) = mb(&k, &an);
            let bs = shear_b(&k, xi, eta, &tied);
            gemm_bt_d_b(&mut ke, nd, &bm, 3, &law.a, k.da);
            gemm_bt_d_b(&mut ke, nd, &bb, 3, &law.d, k.da);
            gemm_bt_d_b(&mut ke, nd, &bs, 2, &law.ds, k.da);
            area += k.da;
            e3_avg = add(e3_avg, k.e3);
        }
    }
    let e3 = if norm(e3_avg) > 1e-14 {
        normalize(e3_avg)?
    } else {
        an[0]
    };
    let kd = DRILL * k6rot * law.drill_eh * area / 4.0;
    for a in 0..4 {
        for i in 0..3 {
            for j in 0..3 {
                ke[(6 * a + 3 + i) * nd + (6 * a + 3 + j)] += kd * e3[i] * e3[j];
            }
        }
    }
    Ok((ke, area))
}

pub fn s4_stress(
    xyz: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    h: f64,
    directors: Option<&[[f64; 3]]>,
) -> Result<Vec<[f64; 6]>> {
    if xyz.len() < 4 {
        return err("S4: zu wenige Knoten.");
    }
    let mut p = [[0.0; 3]; 4];
    p.copy_from_slice(&xyz[..4]);
    let an = resolve_directors(&p, directors)?;
    let tied = tied_rows(&p, &an)?;
    let dm0 = d_plane_stress(e, nu)?;
    let gk = e / (2.0 * (1.0 + nu)) * K_SHEAR;
    let mut out = vec![[0.0; 6]; 4];
    for a in 0..4 {
        let xi = QUAD_XI[a][0];
        let eta = QUAD_XI[a][1];
        let k = kin(&p, &an, xi, eta)?;
        let (bm, bb) = mb(&k, &an);
        let bs = shear_b(&k, xi, eta, &tied);
        let mut eps = [0.0; 3];
        let mut kap = [0.0; 3];
        let mut gam = [0.0; 2];
        for col in 0..24 {
            let u = ue.get(col).copied().unwrap_or(0.0);
            for i in 0..3 {
                eps[i] += bm[i * 24 + col] * u;
                kap[i] += bb[i * 24 + col] * u;
            }
            gam[0] += bs[col] * u;
            gam[1] += bs[24 + col] * u;
        }
        let z = 0.5 * h;
        let em = [eps[0] + z * kap[0], eps[1] + z * kap[1], eps[2] + z * kap[2]];
        let mut sm = [0.0; 3];
        for i in 0..3 {
            for j in 0..3 {
                sm[i] += dm0[i * 3 + j] * em[j];
            }
        }
        let sl = [sm[0], sm[1], 0.0, sm[2], gk * gam[1], gk * gam[0]];
        out[a] = tensor_rotate(sl, k.e1, k.e2, k.e3);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flat_pure_bending_energy() {
        let xyz = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let kappa = 0.1;
        let mut ue = vec![0.0; 24];
        for a in 0..4 {
            let x = xyz[a][0];
            ue[6 * a + 2] = -0.5 * kappa * x * x;
            ue[6 * a + 4] = kappa * x;
        }
        let (ke, _) = s4_ke(&xyz, 2.1e11, 0.3, 0.01, None).unwrap();
        let mut e = 0.0;
        for i in 0..24 {
            let mut ku = 0.0;
            for j in 0..24 {
                ku += ke[i * 24 + j] * ue[j];
            }
            e += ue[i] * ku;
        }
        e *= 0.5;
        // D = E h^3 / 12(1-ν²) = 19230.8, ½ D κ² A = 192.31
        assert!((e - 192.307).abs() < 1.0, "Ebnd={e}");
    }

    /// NAFEMS LE3, quarter hemisphere with 18° cut-out. Target ux(A) = 0.185 m.
    /// Bathe MITC4 on the regular mesh: 1.04 / 1.01 / 1.00 at 4×4 / 8×8 / 16×16.
    #[test]
    fn le3_matches_bathe() {
        let target = 0.185;
        for n in [4usize, 8] {
            let (inp, a_id) = le3_cutout(n);
            let out = crate::solve_native(&inp).unwrap();
            let ux = out.u[out.model.node_index(a_id).unwrap()][0];
            let rel = (ux - target).abs() / target;
            assert!(
                rel < if n == 4 { 0.05 } else { 0.02 },
                "LE3 n={n} ux={ux} rel={rel}"
            );
        }
    }

    fn le3_cutout(n: usize) -> (String, i32) {
        let r = 10.0f64;
        let th0 = 18.0_f64.to_radians();
        let th1 = std::f64::consts::FRAC_PI_2;
        let mut s = String::from("*HEADING\nLE3\n*NODE\n");
        let mut nid = vec![vec![0i32; n + 1]; n + 1];
        let mut id = 0i32;
        for i in 0..=n {
            let th = th0 + (th1 - th0) * (i as f64 / n as f64);
            for j in 0..=n {
                let ph = std::f64::consts::FRAC_PI_2 * (j as f64 / n as f64);
                id += 1;
                nid[i][j] = id;
                s.push_str(&format!(
                    "{id}, {}, {}, {}\n",
                    r * th.sin() * ph.cos(),
                    r * th.sin() * ph.sin(),
                    r * th.cos()
                ));
            }
        }
        s.push_str("*ELEMENT, TYPE=S4, ELSET=HEMI\n");
        let mut eid = 0i32;
        for i in 0..n {
            for j in 0..n {
                eid += 1;
                s.push_str(&format!(
                    "{eid}, {}, {}, {}, {}\n",
                    nid[i][j],
                    nid[i + 1][j],
                    nid[i + 1][j + 1],
                    nid[i][j + 1]
                ));
            }
        }
        s.push_str("*NSET, NSET=AE\n");
        for i in 0..=n {
            s.push_str(&format!("{}\n", nid[i][0]));
        }
        s.push_str("*NSET, NSET=CE\n");
        for i in 0..=n {
            s.push_str(&format!("{}\n", nid[i][n]));
        }
        let a = nid[n][0];
        let c = nid[n][n];
        let ept = nid[0][0];
        s.push_str(&format!(
            "*MATERIAL, NAME=AL\n*ELASTIC\n6.825e10, 0.3\n*SHELL SECTION, ELSET=HEMI, MATERIAL=AL\n0.04\n*BOUNDARY\nAE, 2, 2\nAE, 4, 4\nAE, 6, 6\nCE, 1, 1\nCE, 5, 5\nCE, 6, 6\n{ept}, 3, 3\n*STEP\n*STATIC\n*CLOAD\n{a}, 1, 2000\n{c}, 2, -2000\n*END STEP\n"
        ));
        (s, a)
    }
}
