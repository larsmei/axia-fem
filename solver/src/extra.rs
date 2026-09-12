//! Extra continuum / discrete elements: C3D6, C3D8I, C3D8R, T3D2/T3D3, SPRINGA.

use crate::elem::{
    d_iso_3d, fill_b3, gemm_bt_d_b, invert3, sigma_from_b, G2, hex8_dndx, hex8_shape, QUAD_XI,
};
use crate::error::{err, Result};

fn invert_n(a: &[f64], n: usize) -> Result<Vec<f64>> {
    let mut m = a.to_vec();
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        inv[i * n + i] = 1.0;
    }
    for k in 0..n {
        let mut piv = k;
        let mut best = m[k * n + k].abs();
        for i in k + 1..n {
            let v = m[i * n + k].abs();
            if v > best {
                best = v;
                piv = i;
            }
        }
        if best < 1e-18 {
            return err("Kleine Matrix singulär (inkompatible Moden / Kondensation).");
        }
        if piv != k {
            for j in 0..n {
                m.swap(k * n + j, piv * n + j);
                inv.swap(k * n + j, piv * n + j);
            }
        }
        let d = m[k * n + k];
        for j in 0..n {
            m[k * n + j] /= d;
            inv[k * n + j] /= d;
        }
        for i in 0..n {
            if i == k {
                continue;
            }
            let f = m[i * n + k];
            for j in 0..n {
                m[i * n + j] -= f * m[k * n + j];
                inv[i * n + j] -= f * inv[k * n + j];
            }
        }
    }
    Ok(inv)
}

fn hex_as_8(xyz: &[[f64; 3]]) -> [[f64; 3]; 8] {
    let mut p = [[0.0; 3]; 8];
    for i in 0..8 {
        p[i] = xyz[i];
    }
    p
}

// ---------------------------------------------------------------------------
// C3D6 wedge (6-node pentahedron)
// ---------------------------------------------------------------------------

fn wedge_shape(xi: f64, eta: f64, zeta: f64) -> ([f64; 6], [[f64; 3]; 6]) {
    let l1 = xi;
    let l2 = eta;
    let l3 = 1.0 - xi - eta;
    let zm = 0.5 * (1.0 - zeta);
    let zp = 0.5 * (1.0 + zeta);
    let n = [
        l1 * zm,
        l2 * zm,
        l3 * zm,
        l1 * zp,
        l2 * zp,
        l3 * zp,
    ];
    let mut dn = [[0.0; 3]; 6];
    // d/dξ, d/dη, d/dζ
    dn[0] = [zm, 0.0, -0.5 * l1];
    dn[1] = [0.0, zm, -0.5 * l2];
    dn[2] = [-zm, -zm, -0.5 * l3];
    dn[3] = [zp, 0.0, 0.5 * l1];
    dn[4] = [0.0, zp, 0.5 * l2];
    dn[5] = [-zp, -zp, 0.5 * l3];
    (n, dn)
}

pub(crate) fn wedge_dndx(xyz: &[[f64; 3]], xi: f64, eta: f64, zeta: f64) -> Result<([[f64; 3]; 6], f64, [f64; 6])> {
    let (n, dn) = wedge_shape(xi, eta, zeta);
    let mut j = [[0.0; 3]; 3];
    for a in 0..6 {
        for p in 0..3 {
            for q in 0..3 {
                j[q][p] += dn[a][p] * xyz[a][q];
            }
        }
    }
    let (inv, det) = invert3(j)?;
    let mut dndx = [[0.0; 3]; 6];
    for a in 0..6 {
        for i in 0..3 {
            dndx[a][i] = inv[0][i] * dn[a][0] + inv[1][i] * dn[a][1] + inv[2][i] * dn[a][2];
        }
    }
    Ok((dndx, det, n))
}

pub fn wedge6_stiffness(xyz: &[[f64; 3]], e: f64, nu: f64) -> Result<(Vec<f64>, f64)> {
    let d = d_iso_3d(e, nu)?;
    let nd = 18usize;
    let mut ke = vec![0.0; nd * nd];
    let mut vol = 0.0;
    // triangle 3-point (area 1/2) × ζ ±1/√3
    let tri = [[1.0 / 6.0, 1.0 / 6.0], [2.0 / 3.0, 1.0 / 6.0], [1.0 / 6.0, 2.0 / 3.0]];
    let wtri = 1.0 / 6.0;
    let zpts = [-G2, G2];
    for t in &tri {
        for &zeta in &zpts {
            let (dndx, det, _) = wedge_dndx(xyz, t[0], t[1], zeta)?;
            if det <= 0.0 {
                return err("C3D6: negative Jakobideterminante (Knotenreihenfolge).");
            }
            let w = det * wtri;
            let mut b = vec![0.0; 6 * nd];
            fill_b3(&mut b, 6, &dndx);
            gemm_bt_d_b(&mut ke, nd, &b, 6, &d, w);
            vol += w;
        }
    }
    Ok((ke, vol))
}

pub fn wedge6_nodal_stress(xyz: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<Vec<[f64; 6]>> {
    let d = d_iso_3d(e, nu)?;
    let nd = 18usize;
    let corners = [
        [1.0, 0.0, -1.0],
        [0.0, 1.0, -1.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0],
    ];
    let mut out = vec![[0.0; 6]; 6];
    for a in 0..6 {
        let (dndx, _, _) = wedge_dndx(xyz, corners[a][0], corners[a][1], corners[a][2])?;
        let mut b = vec![0.0; 6 * nd];
        fill_b3(&mut b, 6, &dndx);
        let s = sigma_from_b(&b, 6, nd, &d, ue);
        out[a].copy_from_slice(&s);
    }
    Ok(out)
}

pub fn wedge6_body_force(xyz: &[[f64; 3]], bx: f64, by: f64, bz: f64) -> Result<Vec<f64>> {
    let mut fe = vec![0.0; 18];
    let tri = [[1.0 / 6.0, 1.0 / 6.0], [2.0 / 3.0, 1.0 / 6.0], [1.0 / 6.0, 2.0 / 3.0]];
    let wtri = 1.0 / 6.0;
    let zpts = [-G2, G2];
    for t in &tri {
        for &zeta in &zpts {
            let (n, det, _) = {
                let (n, _, det) = {
                    let r = wedge_dndx(xyz, t[0], t[1], zeta)?;
                    (r.2, r.0, r.1)
                };
                (n, det, 0)
            };
            if det <= 0.0 {
                return err("C3D6: negative Jakobideterminante.");
            }
            let w = det * wtri;
            for a in 0..6 {
                fe[3 * a] += n[a] * bx * w;
                fe[3 * a + 1] += n[a] * by * w;
                fe[3 * a + 2] += n[a] * bz * w;
            }
        }
    }
    Ok(fe)
}

pub fn wedge6_face_pressure(xyz: &[[f64; 3]], face: i32, p: f64) -> Result<Vec<f64>> {
    let mut fe = vec![0.0; 18];
    match face {
        1 => tri_face_load(&mut fe, xyz, [0, 1, 2], p),
        2 => tri_face_load(&mut fe, xyz, [3, 5, 4], p), // outward: 4-6-5 → 3,5,4
        3 => quad_face_load(&mut fe, xyz, [0, 1, 4, 3], p)?,
        4 => quad_face_load(&mut fe, xyz, [1, 2, 5, 4], p)?,
        5 => quad_face_load(&mut fe, xyz, [2, 0, 3, 5], p)?,
        _ => return err(format!("Ungültige C3D6-Fläche P{face}")),
    }
    Ok(fe)
}

fn tri_face_load(fe: &mut [f64], xyz: &[[f64; 3]], idx: [usize; 3], p: f64) {
    let a = xyz[idx[0]];
    let b = xyz[idx[1]];
    let c = xyz[idx[2]];
    let ux = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let vx = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let nx = ux[1] * vx[2] - ux[2] * vx[1];
    let ny = ux[2] * vx[0] - ux[0] * vx[2];
    let nz = ux[0] * vx[1] - ux[1] * vx[0];
    // area vector = 0.5 n, consistent load 1/3 each, traction = -p n_hat * area = -p * 0.5 n
    let s = -p * 0.5 / 3.0;
    for &q in &idx {
        fe[3 * q] += s * nx;
        fe[3 * q + 1] += s * ny;
        fe[3 * q + 2] += s * nz;
    }
}

fn quad_face_load(fe: &mut [f64], xyz: &[[f64; 3]], idx: [usize; 4], p: f64) -> Result<()> {
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
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
            let mut rxi = [0.0; 3];
            let mut reta = [0.0; 3];
            for a in 0..4 {
                let q = xyz[idx[a]];
                for k in 0..3 {
                    rxi[k] += dnxi[a] * q[k];
                    reta[k] += dneta[a] * q[k];
                }
            }
            let nx = rxi[1] * reta[2] - rxi[2] * reta[1];
            let ny = rxi[2] * reta[0] - rxi[0] * reta[2];
            let nz = rxi[0] * reta[1] - rxi[1] * reta[0];
            for a in 0..4 {
                let q = idx[a];
                fe[3 * q] += -p * n[a] * nx;
                fe[3 * q + 1] += -p * n[a] * ny;
                fe[3 * q + 2] += -p * n[a] * nz;
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// C3D8I — Wilson/Taylor incompatible modes, 9 internal DOFs condensed
// ---------------------------------------------------------------------------

pub fn hex8i_stiffness(xyz: &[[f64; 3]], e: f64, nu: f64) -> Result<(Vec<f64>, f64)> {
    let p = hex_as_8(xyz);
    let dmat = d_iso_3d(e, nu)?;
    let n = 24usize;
    let ni = 9usize;
    let mut kee = vec![0.0; n * n];
    let mut kea = vec![0.0; n * ni];
    let mut kaa = vec![0.0; ni * ni];
    let mut vol = 0.0;

    let (_, det0, _) = hex8_dndx(&p, 0.0, 0.0, 0.0)?;
    if det0 <= 0.0 {
        return err("C3D8I: negative Jakobideterminante.");
    }
    let (_, dn0) = hex8_shape(0.0, 0.0, 0.0);
    let mut j0 = [[0.0; 3]; 3];
    for a in 0..8 {
        for pidx in 0..3 {
            for q in 0..3 {
                j0[q][pidx] += dn0[a][pidx] * p[a][q];
            }
        }
    }
    let (inv0, _) = invert3(j0)?;

    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            for &zeta in &pts {
                let (dndx, det, _) = hex8_dndx(&p, xi, eta, zeta)?;
                if det <= 0.0 {
                    return err("C3D8I: negative Jakobideterminante.");
                }
                let mut be = vec![0.0; 6 * n];
                fill_b3(&mut be, 8, &dndx);
                gemm_bt_d_b(&mut kee, n, &be, 6, &dmat, det);
                vol += det;

                // incompatible modes mapped with J at center (Taylor)
                let mut ba = vec![0.0; 6 * ni];
                let dparent = [
                    [-2.0 * xi, 0.0, 0.0],
                    [0.0, -2.0 * eta, 0.0],
                    [0.0, 0.0, -2.0 * zeta],
                ];
                for m in 0..3 {
                    let mut g = [0.0; 3];
                    for i in 0..3 {
                        g[i] = inv0[0][i] * dparent[m][0]
                            + inv0[1][i] * dparent[m][1]
                            + inv0[2][i] * dparent[m][2];
                    }
                    for dir in 0..3 {
                        let col = m * 3 + dir;
                        // B column for displacement in `dir`
                        match dir {
                            0 => {
                                ba[0 * ni + col] = g[0];
                                ba[3 * ni + col] = g[1];
                                ba[5 * ni + col] = g[2];
                            }
                            1 => {
                                ba[1 * ni + col] = g[1];
                                ba[3 * ni + col] = g[0];
                                ba[4 * ni + col] = g[2];
                            }
                            _ => {
                                ba[2 * ni + col] = g[2];
                                ba[4 * ni + col] = g[1];
                                ba[5 * ni + col] = g[0];
                            }
                        }
                    }
                }
                // tmp = D * Ba (6 x 9)
                let mut tmp = vec![0.0; 6 * ni];
                for i in 0..6 {
                    for j in 0..ni {
                        let mut s = 0.0;
                        for k in 0..6 {
                            s += dmat[i * 6 + k] * ba[k * ni + j];
                        }
                        tmp[i * ni + j] = s;
                    }
                }
                // kaa += w Ba^T tmp
                for i in 0..ni {
                    for j in 0..ni {
                        let mut s = 0.0;
                        for k in 0..6 {
                            s += ba[k * ni + i] * tmp[k * ni + j];
                        }
                        kaa[i * ni + j] += det * s;
                    }
                }
                // kea += w Be^T tmp   (24 x 9)
                for i in 0..n {
                    for j in 0..ni {
                        let mut s = 0.0;
                        for k in 0..6 {
                            s += be[k * n + i] * tmp[k * ni + j];
                        }
                        kea[i * ni + j] += det * s;
                    }
                }
            }
        }
    }
    let kaa_inv = invert_n(&kaa, ni)?;
    // Kee -= Kea * Kaa^{-1} * Kea^T
    // t = Kaa^{-1} * Kea^T  (9 x 24)
    let mut t = vec![0.0; ni * n];
    for j in 0..n {
        for i in 0..ni {
            let mut s = 0.0;
            for k in 0..ni {
                s += kaa_inv[i * ni + k] * kea[j * ni + k];
            }
            t[i * n + j] = s;
        }
    }
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in 0..ni {
                s += kea[i * ni + k] * t[k * n + j];
            }
            kee[i * n + j] -= s;
        }
    }
    Ok((kee, vol))
}


pub fn hex8i_nodal_stress(xyz: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<Vec<[f64; 6]>> {
    // Compatible B at nodes is a reasonable recovery for incompatible-mode elements.
    let p = hex_as_8(xyz);
    crate::elem::hex8_nodal_stress(&p, ue, e, nu)
}

// ---------------------------------------------------------------------------
// C3D8R — 1-point + Flanagan–Belytschko hourglass
// ---------------------------------------------------------------------------

pub fn hex8r_stiffness(xyz: &[[f64; 3]], e: f64, nu: f64) -> Result<(Vec<f64>, f64)> {
    let p = hex_as_8(xyz);
    let dmat = d_iso_3d(e, nu)?;
    let n = 24usize;
    let mut ke = vec![0.0; n * n];
    let (dndx, det, _) = hex8_dndx(&p, 0.0, 0.0, 0.0)?;
    if det <= 0.0 {
        return err("C3D8R: negative Jakobideterminante.");
    }
    let vol = 8.0 * det;
    let mut b = vec![0.0; 6 * n];
    fill_b3(&mut b, 8, &dndx);
    gemm_bt_d_b(&mut ke, n, &b, 6, &dmat, vol);

    // Hourglass base vectors (ξη, ηζ, ζξ, ξηζ)
    let h = [
        [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0],
        [1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0],
        [-1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0],
    ];
    let mut gamma = [[0.0; 8]; 4];
    for alpha in 0..4 {
        let mut hx = [0.0; 3];
        for a in 0..8 {
            for j in 0..3 {
                hx[j] += h[alpha][a] * p[a][j];
            }
        }
        for a in 0..8 {
            let mut g = h[alpha][a];
            for j in 0..3 {
                g -= hx[j] * dndx[a][j];
            }
            gamma[alpha][a] = g;
        }
    }
    let mu = e / (2.0 * (1.0 + nu));
    let coeff = 0.1 * mu * vol.powf(1.0 / 3.0);
    for alpha in 0..4 {
        let mut g2 = 0.0;
        for a in 0..8 {
            g2 += gamma[alpha][a] * gamma[alpha][a];
        }
        if g2 < 1e-30 {
            continue;
        }
        let q = coeff / g2.max(1e-30);
        for a in 0..8 {
            for bnode in 0..8 {
                let v = q * gamma[alpha][a] * gamma[alpha][bnode];
                for dir in 0..3 {
                    ke[(3 * a + dir) * n + (3 * bnode + dir)] += v;
                }
            }
        }
    }
    Ok((ke, vol))
}

pub fn hex8r_nodal_stress(xyz: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<Vec<[f64; 6]>> {
    let p = hex_as_8(xyz);
    let dmat = d_iso_3d(e, nu)?;
    let (dndx, _, _) = hex8_dndx(&p, 0.0, 0.0, 0.0)?;
    let n = 24usize;
    let mut b = vec![0.0; 6 * n];
    fill_b3(&mut b, 8, &dndx);
    let s = sigma_from_b(&b, 6, n, &dmat, ue);
    let mut six = [0.0; 6];
    six.copy_from_slice(&s);
    Ok(vec![six; 8])
}

// ---------------------------------------------------------------------------
// T3D2 / T3D3 truss
// ---------------------------------------------------------------------------

pub fn truss_stiffness(xyz: &[[f64; 3]], e: f64, area: f64) -> Result<(Vec<f64>, f64)> {
    let nn = xyz.len();
    if nn < 2 {
        return err("Fachwerk ohne zwei Knoten.");
    }
    let (i0, i1) = if nn == 2 { (0, 1) } else { (0, 2) }; // T3D3: ends 1 and 3, mid unused structurally as linear
    let mut d = [
        xyz[i1][0] - xyz[i0][0],
        xyz[i1][1] - xyz[i0][1],
        xyz[i1][2] - xyz[i0][2],
    ];
    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if len < 1e-18 {
        return err("T3D2: Länge null.");
    }
    d[0] /= len;
    d[1] /= len;
    d[2] /= len;
    let k = e * area.max(1e-30) / len;
    let nd = 3 * nn;
    let mut ke = vec![0.0; nd * nd];
    // linear: only end nodes; T3D3 mid gets a small axial interpolation (quadratic bar)
    if nn == 2 {
        scatter_nn(&mut ke, nd, 0, 0, k, d);
        scatter_nn(&mut ke, nd, 1, 1, k, d);
        scatter_nn(&mut ke, nd, 0, 1, -k, d);
        scatter_nn(&mut ke, nd, 1, 0, -k, d);
    } else {
        // 3-node quadratic truss, 2-point Gauss on [-1,1]
        // N1 = -0.5 ξ (1-ξ), N2 = 0.5 ξ (1+ξ), N3 = 1-ξ²  (ccx: nodes 1,2 ends, 3 mid)
        let pts = [-G2, G2];
        for &xi in &pts {
            let nshp = [-0.5 * xi * (1.0 - xi), 0.5 * xi * (1.0 + xi), 1.0 - xi * xi];
            let dn = [-0.5 + xi, 0.5 + xi, -2.0 * xi];
            let mut dxdxi = [0.0; 3];
            for a in 0..3 {
                for kdir in 0..3 {
                    dxdxi[kdir] += dn[a] * xyz[a][kdir];
                }
            }
            let jac = (dxdxi[0] * dxdxi[0] + dxdxi[1] * dxdxi[1] + dxdxi[2] * dxdxi[2]).sqrt();
            if jac < 1e-18 {
                return err("T3D3: singulär.");
            }
            let mut tdir = [0.0; 3];
            for kdir in 0..3 {
                tdir[kdir] = dxdxi[kdir] / jac;
            }
            let mut dndx_ax = [0.0; 3];
            for a in 0..3 {
                dndx_ax[a] = dn[a] / jac;
            }
            let ea = e * area * jac; // weight = jac * 1
            for a in 0..3 {
                for bnode in 0..3 {
                    let c = ea * dndx_ax[a] * dndx_ax[bnode];
                    for i in 0..3 {
                        for j in 0..3 {
                            ke[(3 * a + i) * nd + (3 * bnode + j)] += c * tdir[i] * tdir[j];
                        }
                    }
                }
            }
            let _ = nshp;
        }
    }
    Ok((ke, len * area))
}

fn scatter_nn(ke: &mut [f64], nd: usize, a: usize, b: usize, k: f64, n: [f64; 3]) {
    for i in 0..3 {
        for j in 0..3 {
            ke[(3 * a + i) * nd + (3 * b + j)] += k * n[i] * n[j];
        }
    }
}

pub fn truss_nodal_stress(xyz: &[[f64; 3]], ue: &[f64], e: f64, area: f64) -> Result<Vec<[f64; 6]>> {
    let nn = xyz.len();
    let i1 = if nn == 2 { 1 } else { 2 };
    let mut d = [
        xyz[i1][0] - xyz[0][0],
        xyz[i1][1] - xyz[0][1],
        xyz[i1][2] - xyz[0][2],
    ];
    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt().max(1e-18);
    d[0] /= len;
    d[1] /= len;
    d[2] /= len;
    let du = [
        ue[3 * i1] - ue[0],
        ue[3 * i1 + 1] - ue[1],
        ue[3 * i1 + 2] - ue[2],
    ];
    let eps = (du[0] * d[0] + du[1] * d[1] + du[2] * d[2]) / len;
    let sig = e * eps;
    let s = [
        sig * d[0] * d[0],
        sig * d[1] * d[1],
        sig * d[2] * d[2],
        sig * d[0] * d[1],
        sig * d[1] * d[2],
        sig * d[2] * d[0],
    ];
    let _ = area;
    Ok(vec![s; nn])
}

pub fn truss_body_force(xyz: &[[f64; 3]], area: f64, bx: f64, by: f64, bz: f64) -> Vec<f64> {
    let nn = xyz.len();
    let i1 = if nn == 2 { 1 } else { 2 };
    let dx = xyz[i1][0] - xyz[0][0];
    let dy = xyz[i1][1] - xyz[0][1];
    let dz = xyz[i1][2] - xyz[0][2];
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    let w = area * len / nn as f64;
    let mut fe = vec![0.0; 3 * nn];
    for a in 0..nn {
        fe[3 * a] = w * bx;
        fe[3 * a + 1] = w * by;
        fe[3 * a + 2] = w * bz;
    }
    fe
}

// ---------------------------------------------------------------------------
// SPRINGA — axial spring between two nodes
// ---------------------------------------------------------------------------

pub fn spring_stiffness(xyz: &[[f64; 3]], k: f64) -> Result<(Vec<f64>, f64)> {
    let mut d = [
        xyz[1][0] - xyz[0][0],
        xyz[1][1] - xyz[0][1],
        xyz[1][2] - xyz[0][2],
    ];
    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if len < 1e-18 {
        d = [1.0, 0.0, 0.0];
    } else {
        d[0] /= len;
        d[1] /= len;
        d[2] /= len;
    }
    let mut ke = vec![0.0; 36];
    scatter_nn(&mut ke, 6, 0, 0, k, d);
    scatter_nn(&mut ke, 6, 1, 1, k, d);
    scatter_nn(&mut ke, 6, 0, 1, -k, d);
    scatter_nn(&mut ke, 6, 1, 0, -k, d);
    Ok((ke, len))
}

pub fn truss_thermal_force(xyz: &[[f64; 3]], e: f64, area: f64, alpha: f64, dt: f64) -> Vec<f64> {
    let nn = xyz.len();
    let i1 = if nn == 2 { 1 } else { nn - 1 };
    let mut d = [
        xyz[i1][0] - xyz[0][0],
        xyz[i1][1] - xyz[0][1],
        xyz[i1][2] - xyz[0][2],
    ];
    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt().max(1e-18);
    d[0] /= len;
    d[1] /= len;
    d[2] /= len;
    let n = e * area * alpha * dt;
    let mut fe = vec![0.0; 3 * nn];
    for k in 0..3 {
        fe[k] -= n * d[k];
        fe[3 * i1 + k] += n * d[k];
    }
    fe
}

pub fn hex8_thermal_force(xyz: &[[f64; 3]], e: f64, nu: f64, alpha: f64, dt: f64) -> Result<Vec<f64>> {
    let p = hex_as_8(xyz);
    let dmat = d_iso_3d(e, nu)?;
    let mut eps = [0.0; 6];
    eps[0] = alpha * dt;
    eps[1] = alpha * dt;
    eps[2] = alpha * dt;
    let mut sig = [0.0; 6];
    for i in 0..6 {
        for j in 0..6 {
            sig[i] += dmat[i * 6 + j] * eps[j];
        }
    }
    let n = 24usize;
    let mut fe = vec![0.0; n];
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            for &zeta in &pts {
                let (dndx, det, _) = hex8_dndx(&p, xi, eta, zeta)?;
                if det <= 0.0 {
                    continue;
                }
                let mut b = vec![0.0; 6 * n];
                fill_b3(&mut b, 8, &dndx);
                for j in 0..n {
                    let mut s = 0.0;
                    for i in 0..6 {
                        s += b[i * n + j] * sig[i];
                    }
                    fe[j] += s * det;
                }
            }
        }
    }
    Ok(fe)
}

pub fn spring_nodal_stress(_xyz: &[[f64; 3]], _ue: &[f64], _k: f64) -> Vec<[f64; 6]> {
    vec![[0.0; 6]; 2]
}

// ---------------------------------------------------------------------------
// C3D15 — 15-node quadratic wedge (pentahedron)
// Nodes: 1-3 bottom corners, 4-6 top corners,
// 7-9 bottom mids (1-2,2-3,3-1), 10-12 top mids, 13-15 vertical mids.
// Parent: L1=ξ, L2=η, L3=1-ξ-η, ζ ∈ [-1, 1].
// ---------------------------------------------------------------------------

const G3: [f64; 3] = [-0.7745966692414834, 0.0, 0.7745966692414834];
const W3: [f64; 3] = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];

fn wedge15_shape(xi: f64, eta: f64, zeta: f64) -> ([f64; 15], [[f64; 3]; 15]) {
    let l1 = xi;
    let l2 = eta;
    let l3 = 1.0 - xi - eta;
    let zm = 1.0 - zeta;
    let zp = 1.0 + zeta;
    let mut n = [0.0; 15];
    let mut dn = [[0.0; 3]; 15];

    // corners bottom 1-3, top 4-6
    n[0] = 0.5 * l1 * zm * (2.0 * l1 - 2.0 - zeta);
    n[1] = 0.5 * l2 * zm * (2.0 * l2 - 2.0 - zeta);
    n[2] = 0.5 * l3 * zm * (2.0 * l3 - 2.0 - zeta);
    n[3] = 0.5 * l1 * zp * (2.0 * l1 - 2.0 + zeta);
    n[4] = 0.5 * l2 * zp * (2.0 * l2 - 2.0 + zeta);
    n[5] = 0.5 * l3 * zp * (2.0 * l3 - 2.0 + zeta);
    // bottom mids 7-9, top mids 10-12
    n[6] = 2.0 * l1 * l2 * zm;
    n[7] = 2.0 * l2 * l3 * zm;
    n[8] = 2.0 * l3 * l1 * zm;
    n[9] = 2.0 * l1 * l2 * zp;
    n[10] = 2.0 * l2 * l3 * zp;
    n[11] = 2.0 * l3 * l1 * zp;
    // vertical mids 13-15
    n[12] = l1 * (1.0 - zeta * zeta);
    n[13] = l2 * (1.0 - zeta * zeta);
    n[14] = l3 * (1.0 - zeta * zeta);

    let d1b = 0.5 * zm * (4.0 * l1 - 2.0 - zeta);
    let d2b = 0.5 * zm * (4.0 * l2 - 2.0 - zeta);
    let d3b = 0.5 * zm * (4.0 * l3 - 2.0 - zeta);
    let d1t = 0.5 * zp * (4.0 * l1 - 2.0 + zeta);
    let d2t = 0.5 * zp * (4.0 * l2 - 2.0 + zeta);
    let d3t = 0.5 * zp * (4.0 * l3 - 2.0 + zeta);

    // d/dξ, d/dη, d/dζ
    dn[0] = [d1b, 0.0, -0.5 * l1 * (2.0 * l1 - 1.0 - 2.0 * zeta)];
    dn[1] = [0.0, d2b, -0.5 * l2 * (2.0 * l2 - 1.0 - 2.0 * zeta)];
    dn[2] = [-d3b, -d3b, -0.5 * l3 * (2.0 * l3 - 1.0 - 2.0 * zeta)];
    dn[3] = [d1t, 0.0, 0.5 * l1 * (2.0 * l1 - 1.0 + 2.0 * zeta)];
    dn[4] = [0.0, d2t, 0.5 * l2 * (2.0 * l2 - 1.0 + 2.0 * zeta)];
    dn[5] = [-d3t, -d3t, 0.5 * l3 * (2.0 * l3 - 1.0 + 2.0 * zeta)];

    dn[6] = [2.0 * l2 * zm, 2.0 * l1 * zm, -2.0 * l1 * l2];
    dn[7] = [-2.0 * l2 * zm, 2.0 * (l3 - l2) * zm, -2.0 * l2 * l3];
    dn[8] = [2.0 * (l3 - l1) * zm, -2.0 * l1 * zm, -2.0 * l3 * l1];
    dn[9] = [2.0 * l2 * zp, 2.0 * l1 * zp, 2.0 * l1 * l2];
    dn[10] = [-2.0 * l2 * zp, 2.0 * (l3 - l2) * zp, 2.0 * l2 * l3];
    dn[11] = [2.0 * (l3 - l1) * zp, -2.0 * l1 * zp, 2.0 * l3 * l1];

    let o = 1.0 - zeta * zeta;
    dn[12] = [o, 0.0, -2.0 * zeta * l1];
    dn[13] = [0.0, o, -2.0 * zeta * l2];
    dn[14] = [-o, -o, -2.0 * zeta * l3];
    (n, dn)
}

pub(crate) fn wedge15_dndx(
    xyz: &[[f64; 3]],
    xi: f64,
    eta: f64,
    zeta: f64,
) -> Result<([[f64; 3]; 15], f64, [f64; 15])> {
    let (n, dn) = wedge15_shape(xi, eta, zeta);
    let mut j = [[0.0; 3]; 3];
    for a in 0..15 {
        for p in 0..3 {
            for q in 0..3 {
                j[q][p] += dn[a][p] * xyz[a][q];
            }
        }
    }
    let (inv, det) = invert3(j)?;
    let mut dndx = [[0.0; 3]; 15];
    for a in 0..15 {
        for i in 0..3 {
            dndx[a][i] = inv[0][i] * dn[a][0] + inv[1][i] * dn[a][1] + inv[2][i] * dn[a][2];
        }
    }
    Ok((dndx, det, n))
}

fn wedge15_gauss() -> Vec<(f64, f64, f64, f64)> {
    let tri = [
        [1.0 / 6.0, 1.0 / 6.0],
        [2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0],
    ];
    let wtri = 1.0 / 6.0;
    let mut o = Vec::new();
    for t in &tri {
        for k in 0..3 {
            o.push((t[0], t[1], G3[k], wtri * W3[k]));
        }
    }
    o
}

pub fn wedge15_stiffness(xyz: &[[f64; 3]], e: f64, nu: f64) -> Result<(Vec<f64>, f64)> {
    if xyz.len() < 15 {
        return err("C3D15 braucht 15 Knoten.");
    }
    let d = d_iso_3d(e, nu)?;
    let nd = 45usize;
    let mut ke = vec![0.0; nd * nd];
    let mut vol = 0.0;
    for (xi, eta, zeta, w0) in wedge15_gauss() {
        let (dndx, det, _) = wedge15_dndx(xyz, xi, eta, zeta)?;
        if det <= 0.0 {
            return err("C3D15: negative Jakobideterminante (Knotenreihenfolge).");
        }
        let w = det * w0;
        let mut b = vec![0.0; 6 * nd];
        fill_b3(&mut b, 15, &dndx);
        gemm_bt_d_b(&mut ke, nd, &b, 6, &d, w);
        vol += w;
    }
    Ok((ke, vol))
}

pub fn wedge15_nodal_stress(xyz: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<Vec<[f64; 6]>> {
    let d = d_iso_3d(e, nu)?;
    let nd = 45usize;
    let corners = [
        [1.0, 0.0, -1.0],
        [0.0, 1.0, -1.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0],
        [0.5, 0.5, -1.0],
        [0.0, 0.5, -1.0],
        [0.5, 0.0, -1.0],
        [0.5, 0.5, 1.0],
        [0.0, 0.5, 1.0],
        [0.5, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
    ];
    let mut out = vec![[0.0; 6]; 15];
    for a in 0..15 {
        let (dndx, _, _) = wedge15_dndx(xyz, corners[a][0], corners[a][1], corners[a][2])?;
        let mut b = vec![0.0; 6 * nd];
        fill_b3(&mut b, 15, &dndx);
        let s = sigma_from_b(&b, 6, nd, &d, ue);
        out[a].copy_from_slice(&s);
    }
    Ok(out)
}

pub fn wedge15_body_force(xyz: &[[f64; 3]], bx: f64, by: f64, bz: f64) -> Result<Vec<f64>> {
    let mut fe = vec![0.0; 45];
    for (xi, eta, zeta, w0) in wedge15_gauss() {
        let (_, det, n) = wedge15_dndx(xyz, xi, eta, zeta)?;
        if det <= 0.0 {
            continue;
        }
        let w = det * w0;
        for a in 0..15 {
            fe[3 * a] += n[a] * bx * w;
            fe[3 * a + 1] += n[a] * by * w;
            fe[3 * a + 2] += n[a] * bz * w;
        }
    }
    Ok(fe)
}

pub fn wedge15_face_pressure(xyz: &[[f64; 3]], face: i32, p: f64) -> Result<Vec<f64>> {
    let mut fe = vec![0.0; 45];
    match face {
        1 => tri6_face_load(&mut fe, xyz, [0, 1, 2, 6, 7, 8], p)?,
        2 => tri6_face_load(&mut fe, xyz, [3, 5, 4, 11, 10, 9], p)?,
        3 => quad8_face_load(&mut fe, xyz, [0, 1, 4, 3, 6, 13, 9, 12], p)?,
        4 => quad8_face_load(&mut fe, xyz, [1, 2, 5, 4, 7, 14, 10, 13], p)?,
        5 => quad8_face_load(&mut fe, xyz, [2, 0, 3, 5, 8, 12, 11, 14], p)?,
        _ => return err(format!("Ungültige C3D15-Fläche P{face}")),
    }
    Ok(fe)
}

fn tri6_face_load(fe: &mut [f64], xyz: &[[f64; 3]], idx: [usize; 6], p: f64) -> Result<()> {
    // 3-point triangle, parent (L1,L2), area 1/2
    let gps = [
        [1.0 / 6.0, 1.0 / 6.0],
        [2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0],
    ];
    let w0 = 1.0 / 6.0;
    for g in &gps {
        let l1 = g[0];
        let l2 = g[1];
        let l3 = 1.0 - l1 - l2;
        let nshp = [
            l1 * (2.0 * l1 - 1.0),
            l2 * (2.0 * l2 - 1.0),
            l3 * (2.0 * l3 - 1.0),
            4.0 * l1 * l2,
            4.0 * l2 * l3,
            4.0 * l3 * l1,
        ];
        let dn1 = [4.0 * l1 - 1.0, 0.0, -(4.0 * l3 - 1.0), 4.0 * l2, -4.0 * l2, 4.0 * (l3 - l1)];
        let dn2 = [0.0, 4.0 * l2 - 1.0, -(4.0 * l3 - 1.0), 4.0 * l1, 4.0 * (l3 - l2), -4.0 * l1];
        let mut rxi = [0.0; 3];
        let mut reta = [0.0; 3];
        for a in 0..6 {
            let q = xyz[idx[a]];
            for k in 0..3 {
                rxi[k] += dn1[a] * q[k];
                reta[k] += dn2[a] * q[k];
            }
        }
        let nx = rxi[1] * reta[2] - rxi[2] * reta[1];
        let ny = rxi[2] * reta[0] - rxi[0] * reta[2];
        let nz = rxi[0] * reta[1] - rxi[1] * reta[0];
        for a in 0..6 {
            let q = idx[a];
            fe[3 * q] += -p * nshp[a] * nx * w0;
            fe[3 * q + 1] += -p * nshp[a] * ny * w0;
            fe[3 * q + 2] += -p * nshp[a] * nz * w0;
        }
    }
    Ok(())
}

fn quad8_face_load(fe: &mut [f64], xyz: &[[f64; 3]], idx: [usize; 8], p: f64) -> Result<()> {
    for i in 0..3 {
        for j in 0..3 {
            let xi = G3[i];
            let eta = G3[j];
            let w = W3[i] * W3[j];
            let (nshp, dn) = crate::quadratic::quad8_shape(xi, eta);
            let mut rxi = [0.0; 3];
            let mut reta = [0.0; 3];
            for a in 0..8 {
                let q = xyz[idx[a]];
                for k in 0..3 {
                    rxi[k] += dn[a][0] * q[k];
                    reta[k] += dn[a][1] * q[k];
                }
            }
            let nx = rxi[1] * reta[2] - rxi[2] * reta[1];
            let ny = rxi[2] * reta[0] - rxi[0] * reta[2];
            let nz = rxi[0] * reta[1] - rxi[1] * reta[0];
            for a in 0..8 {
                let q = idx[a];
                fe[3 * q] += -p * nshp[a] * nx * w;
                fe[3 * q + 1] += -p * nshp[a] * ny * w;
                fe[3 * q + 2] += -p * nshp[a] * nz * w;
            }
        }
    }
    Ok(())
}

pub fn wedge15_kg(xyz: &[[f64; 3]], stress: &[f64; 6]) -> Result<Vec<f64>> {
    let nd = 45usize;
    let mut kg = vec![0.0; nd * nd];
    for (xi, eta, zeta, w0) in wedge15_gauss() {
        let (dndx, det, _) = wedge15_dndx(xyz, xi, eta, zeta)?;
        if det <= 0.0 {
            continue;
        }
        crate::eigen::add_continuum_kg(&mut kg, 15, &dndx, stress, det * w0);
    }
    Ok(kg)
}

pub fn wedge6_kg(xyz: &[[f64; 3]], stress: &[f64; 6]) -> Result<Vec<f64>> {
    let nd = 18usize;
    let mut kg = vec![0.0; nd * nd];
    let tri = [
        [1.0 / 6.0, 1.0 / 6.0],
        [2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0],
    ];
    let wtri = 1.0 / 6.0;
    let zpts = [-G2, G2];
    for t in &tri {
        for &zeta in &zpts {
            let (dndx, det, _) = wedge_dndx(xyz, t[0], t[1], zeta)?;
            if det <= 0.0 {
                continue;
            }
            crate::eigen::add_continuum_kg(&mut kg, 6, &dndx, stress, det * wtri);
        }
    }
    Ok(kg)
}
