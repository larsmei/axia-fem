//! Generalized eigen: *FREQUENCY (K φ = ω² M φ) and *BUCKLE (K φ = λ (−Kg) φ).

use crate::error::{err, Result};
use crate::linalg::{csr_from_triplets, solve_kff};

pub struct EigenResult {
    pub values: Vec<f64>,
    pub vectors: Vec<Vec<f64>>,
}

/// Lowest `nmodes` of K φ = λ M φ with diagonal M, using inverse subspace iteration.
pub fn subspace_gen(
    n: usize,
    k_trips: Vec<(usize, usize, f64)>,
    m_diag: &[f64],
    nmodes: usize,
) -> Result<EigenResult> {
    if n == 0 || nmodes == 0 {
        return err("Eigenproblem ohne freie DOF / ohne Moden.");
    }
    let nvec = nmodes.min(n).saturating_add(4).min(n).max(nmodes.min(n));
    let kcsr = csr_from_triplets(n, k_trips.clone());

    let mut rng = 1u64;
    let mut x = vec![vec![0.0; n]; nvec];
    for v in &mut x {
        for i in 0..n {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            v[i] = ((rng >> 33) as f64 / (1u64 << 31) as f64) - 1.0;
            if m_diag[i] <= 0.0 {
                v[i] = 0.0;
            }
        }
    }
    m_orthonormalize(&mut x, m_diag);

    let mut values = vec![0.0; nvec];
    for _it in 0..40 {
        let mut y = Vec::with_capacity(nvec);
        for v in &x {
            let mut rhs = vec![0.0; n];
            for i in 0..n {
                rhs[i] = m_diag[i] * v[i];
            }
            let s = solve_kff(n, k_trips.clone(), &rhs)?;
            y.push(s.x);
        }
        m_orthonormalize(&mut y, m_diag);
        // Project Kr = Y^T K Y
        let mut kr = vec![0.0; nvec * nvec];
        for j in 0..nvec {
            let mut ky = vec![0.0; n];
            kcsr.matvec(&y[j], &mut ky);
            for i in 0..nvec {
                let mut s = 0.0;
                for k in 0..n {
                    s += y[i][k] * ky[k];
                }
                kr[i * nvec + j] = s;
            }
        }
        let (evals, evecs) = jacobi_eigen(&mut kr, nvec);
        // X = Y * evecs
        let mut nx = vec![vec![0.0; n]; nvec];
        for j in 0..nvec {
            for k in 0..nvec {
                let c = evecs[k * nvec + j];
                if c.abs() == 0.0 {
                    continue;
                }
                for i in 0..n {
                    nx[j][i] += y[k][i] * c;
                }
            }
        }
        m_orthonormalize(&mut nx, m_diag);
        let mut conv = true;
        for i in 0..nvec {
            if (evals[i] - values[i]).abs() > 1e-8 * (1.0 + evals[i].abs()) {
                conv = false;
            }
            values[i] = evals[i];
        }
        x = nx;
        if conv && _it > 4 {
            break;
        }
    }
    // sort ascending
    let mut order: Vec<usize> = (0..nvec).collect();
    order.sort_by(|&a, &b| values[a].partial_cmp(&values[b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut vals = Vec::new();
    let mut vecs = Vec::new();
    for &i in order.iter().take(nmodes.min(nvec)) {
        if values[i].is_finite() && values[i] > 0.0 {
            vals.push(values[i]);
            vecs.push(x[i].clone());
        } else if values[i].is_finite() {
            vals.push(values[i]);
            vecs.push(x[i].clone());
        }
    }
    if vals.is_empty() {
        return err("Keine Eigenwerte gefunden.");
    }
    Ok(EigenResult {
        values: vals,
        vectors: vecs,
    })
}

/// Inverse subspace iteration for K φ = λ A φ (A typically −Kg).
pub fn subspace_ab(
    n: usize,
    k_trips: Vec<(usize, usize, f64)>,
    a_trips: Vec<(usize, usize, f64)>,
    nmodes: usize,
) -> Result<EigenResult> {
    if n == 0 || nmodes == 0 {
        return err("Eigenproblem ohne freie DOF.");
    }
    let nvec = nmodes.min(n).saturating_add(4).min(n).max(nmodes.min(n));
    let acsr = csr_from_triplets(n, a_trips);
    let kcsr = csr_from_triplets(n, k_trips.clone());
    let mut rng = 7u64;
    let mut x = vec![vec![0.0; n]; nvec];
    for v in &mut x {
        for i in 0..n {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            v[i] = ((rng >> 33) as f64 / (1u64 << 31) as f64) - 1.0;
        }
        let nrm = v.iter().map(|z| z * z).sum::<f64>().sqrt().max(1e-30);
        for z in v.iter_mut() {
            *z /= nrm;
        }
    }
    for _ in 0..30 {
        let mut y = Vec::new();
        for v in &x {
            let mut av = vec![0.0; n];
            acsr.matvec(v, &mut av);
            let s = solve_kff(n, k_trips.clone(), &av)?;
            y.push(s.x);
        }
        // Euclidean QR-ish
        for a in 0..y.len() {
            for b in 0..a {
                let dot: f64 = y[a].iter().zip(y[b].iter()).map(|(p, q)| p * q).sum();
                for i in 0..n {
                    y[a][i] -= dot * y[b][i];
                }
            }
            let nrm = y[a].iter().map(|z| z * z).sum::<f64>().sqrt().max(1e-30);
            for z in y[a].iter_mut() {
                *z /= nrm;
            }
        }
        x = y;
    }
    let mut pairs: Vec<(f64, Vec<f64>)> = Vec::new();
    for v in x {
        let mut kv = vec![0.0; n];
        let mut av = vec![0.0; n];
        kcsr.matvec(&v, &mut kv);
        acsr.matvec(&v, &mut av);
        let num: f64 = v.iter().zip(kv.iter()).map(|(p, q)| p * q).sum();
        let den: f64 = v.iter().zip(av.iter()).map(|(p, q)| p * q).sum();
        if den.abs() < 1e-30 {
            continue;
        }
        pairs.push((num / den, v));
    }
    pairs.sort_by(|a, b| a.0.abs().partial_cmp(&b.0.abs()).unwrap_or(std::cmp::Ordering::Equal));
    let take = nmodes.min(pairs.len());
    Ok(EigenResult {
        values: pairs.iter().take(take).map(|p| p.0).collect(),
        vectors: pairs.iter().take(take).map(|p| p.1.clone()).collect(),
    })
}

fn m_orthonormalize(x: &mut [Vec<f64>], m: &[f64]) {
    let n = m.len();
    for a in 0..x.len() {
        for b in 0..a {
            let mut dot = 0.0;
            for i in 0..n {
                dot += x[a][i] * m[i] * x[b][i];
            }
            for i in 0..n {
                x[a][i] -= dot * x[b][i];
            }
        }
        let mut nrm = 0.0;
        for i in 0..n {
            nrm += x[a][i] * m[i] * x[a][i];
        }
        let s = nrm.sqrt();
        if s > 1e-30 {
            for i in 0..n {
                x[a][i] /= s;
            }
        }
    }
}

/// Jacobi eigen-decomposition of a dense symmetric n×n matrix `a` (row-major).
/// Returns (eigenvalues, eigenvectors as columns stored row-major: evecs[k*n + j] is component k of vector j).
fn jacobi_eigen(a: &mut [f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _ in 0..40 * n.max(1) {
        let mut p = 0usize;
        let mut q = 1usize;
        let mut best = 0.0;
        for i in 0..n {
            for j in i + 1..n {
                let x = a[i * n + j].abs();
                if x > best {
                    best = x;
                    p = i;
                    q = j;
                }
            }
        }
        if best < 1e-14 {
            break;
        }
        let app = a[p * n + p];
        let aqq = a[q * n + q];
        let apq = a[p * n + q];
        let tau = (aqq - app) / (2.0 * apq);
        let t = {
            let s = if tau >= 0.0 { 1.0 } else { -1.0 };
            s / (tau.abs() + (1.0 + tau * tau).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        for i in 0..n {
            if i != p && i != q {
                let aip = a[i * n + p];
                let aiq = a[i * n + q];
                a[i * n + p] = c * aip - s * aiq;
                a[p * n + i] = a[i * n + p];
                a[i * n + q] = s * aip + c * aiq;
                a[q * n + i] = a[i * n + q];
            }
        }
        a[p * n + p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[q * n + q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[p * n + q] = 0.0;
        a[q * n + p] = 0.0;
        for i in 0..n {
            let vip = v[i * n + p];
            let viq = v[i * n + q];
            v[i * n + p] = c * vip - s * viq;
            v[i * n + q] = s * vip + c * viq;
        }
    }
    let mut evals = vec![0.0; n];
    for i in 0..n {
        evals[i] = a[i * n + i];
    }
    (evals, v)
}

/// Truss geometric stiffness. `n_axial` tension positive.
pub fn truss_kg(xyz: &[[f64; 3]], n_axial: f64) -> Vec<f64> {
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
    let nd = 3 * nn;
    let mut kg = vec![0.0; nd * nd];
    let c = n_axial / len;
    for a in 0..2 {
        let ia = if a == 0 { 0 } else { i1 };
        for b in 0..2 {
            let ib = if b == 0 { 0 } else { i1 };
            let s = if a == b { c } else { -c };
            for i in 0..3 {
                for j in 0..3 {
                    let proj = if i == j { 1.0 } else { 0.0 } - d[i] * d[j];
                    kg[(3 * ia + i) * nd + (3 * ib + j)] += s * proj;
                }
            }
        }
    }
    kg
}

/// Continuum geometric stiffness increment at one Gauss point:
/// Kg^{ab}_{ii} += (∇N_a · σ · ∇N_b) w  (same for each translation dir).
pub fn add_continuum_kg(kg: &mut [f64], nnode: usize, dndx: &[[f64; 3]], stress: &[f64; 6], w: f64) {
    let n = 3 * nnode;
    let sxx = stress[0];
    let syy = stress[1];
    let szz = stress[2];
    let sxy = stress[3];
    let syz = stress[4];
    let szx = stress[5];
    for a in 0..nnode {
        for b in 0..nnode {
            let ga = dndx[a];
            let gb = dndx[b];
            let gtg = ga[0] * (sxx * gb[0] + sxy * gb[1] + szx * gb[2])
                + ga[1] * (sxy * gb[0] + syy * gb[1] + syz * gb[2])
                + ga[2] * (szx * gb[0] + syz * gb[1] + szz * gb[2]);
            let v = gtg * w;
            for dir in 0..3 {
                kg[(3 * a + dir) * n + (3 * b + dir)] += v;
            }
        }
    }
}

/// Hex8 geometric stiffness from constant stress (centroid).
pub fn hex8_kg(xyz: &[[f64; 3]], stress: &[f64; 6]) -> Result<Vec<f64>> {
    use crate::elem::{hex8_dndx, G2};
    let mut p = [[0.0; 3]; 8];
    for i in 0..8 {
        p[i] = xyz[i];
    }
    let n = 24usize;
    let mut kg = vec![0.0; n * n];
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            for &zeta in &pts {
                let (dndx, det, _) = hex8_dndx(&p, xi, eta, zeta)?;
                if det <= 0.0 {
                    continue;
                }
                add_continuum_kg(&mut kg, 8, &dndx, stress, det);
            }
        }
    }
    Ok(kg)
}

pub fn hex20_kg(xyz: &[[f64; 3]], stress: &[f64; 6], reduced: bool) -> Result<Vec<f64>> {
    let n = 60usize;
    let mut kg = vec![0.0; n * n];
    for (xi, eta, zeta, w) in crate::quadratic::hex_gauss(reduced) {
        let (dndx, det, _) = crate::quadratic::hex20_dndx(xyz, xi, eta, zeta)?;
        if det <= 0.0 {
            continue;
        }
        add_continuum_kg(&mut kg, 20, &dndx, stress, w * det);
    }
    Ok(kg)
}

pub fn tet4_kg(xyz: &[[f64; 3]], stress: &[f64; 6]) -> Result<Vec<f64>> {
    let mut j = [[0.0; 3]; 3];
    for p in 0..3 {
        for q in 0..3 {
            j[q][p] = xyz[p + 1][q] - xyz[0][q];
        }
    }
    let (inv, det) = crate::elem::invert3(j)?;
    if det <= 0.0 {
        return err("C3D4 Kg: negative Jakobideterminante.");
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
    let mut kg = vec![0.0; 12 * 12];
    add_continuum_kg(&mut kg, 4, &dndx, stress, vol);
    Ok(kg)
}

pub fn tet10_kg(xyz: &[[f64; 3]], stress: &[f64; 6]) -> Result<Vec<f64>> {
    let n = 30usize;
    let mut kg = vec![0.0; n * n];
    let a = 0.5854101966249685;
    let b = 0.1381966011250105;
    let w = 1.0 / 24.0;
    let pts = [[b, b, b], [a, b, b], [b, a, b], [b, b, a]];
    for p in &pts {
        let (dndx, det, _) = crate::quadratic::tet10_dndx(xyz, p[0], p[1], p[2])?;
        if det <= 0.0 {
            continue;
        }
        add_continuum_kg(&mut kg, 10, &dndx, stress, w * det);
    }
    Ok(kg)
}

/// Beam geometric stiffness (string-like + bending 2D formula in two planes).
pub fn beam_kg(xyz: &[[f64; 3]], n_axial: f64, nn: usize) -> Result<Vec<f64>> {
    use crate::beam::orthonormal;
    let dx = xyz[xyz.len() - 1][0] - xyz[0][0];
    let dy = xyz[xyz.len() - 1][1] - xyz[0][1];
    let dz = xyz[xyz.len() - 1][2] - xyz[0][2];
    // use first and last of the provided xyz (already only end nodes typically)
    let mut tdir = [xyz[1][0] - xyz[0][0], xyz[1][1] - xyz[0][1], xyz[1][2] - xyz[0][2]];
    if xyz.len() >= 3 {
        tdir = [
            xyz[xyz.len() - 1][0] - xyz[0][0],
            xyz[xyz.len() - 1][1] - xyz[0][1],
            xyz[xyz.len() - 1][2] - xyz[0][2],
        ];
    }
    let len = (tdir[0] * tdir[0] + tdir[1] * tdir[1] + tdir[2] * tdir[2]).sqrt();
    if len < 1e-18 {
        return err("Balken Kg: Länge null");
    }
    tdir[0] /= len;
    tdir[1] /= len;
    tdir[2] /= len;
    let (n1, n2) = orthonormal(tdir, [0.0, 0.0, -1.0])?;
    let nd = 6 * nn;
    let mut kg = vec![0.0; nd * nd];
    let p = n_axial;
    // Local 12x12 for two end nodes, then scatter to ends 0 and last.
    let l = len;
    let c = p / l;
    // translational string on n1 and n2 at both ends
    let ends = [0usize, nn - 1];
    for (ia, &a) in ends.iter().enumerate() {
        for (ib, &b) in ends.iter().enumerate() {
            let s = if ia == ib { c } else { -c };
            for k in 0..3 {
                for m in 0..3 {
                    let proj = n1[k] * n1[m] + n2[k] * n2[m];
                    kg[(6 * a + k) * nd + (6 * b + m)] += s * proj;
                }
            }
        }
    }
    // rotational geometric (2 L/15 style, small)
    let cm = p * l / 15.0;
    for (ia, &a) in ends.iter().enumerate() {
        for (ib, &b) in ends.iter().enumerate() {
            let s = if ia == ib { 2.0 * cm } else { -cm };
            for k in 0..3 {
                for m in 0..3 {
                    let proj = n1[k] * n1[m] + n2[k] * n2[m];
                    kg[(6 * a + 3 + k) * nd + (6 * b + 3 + m)] += s * proj;
                }
            }
        }
    }
    let _ = (dx, dy, dz);
    Ok(kg)
}
