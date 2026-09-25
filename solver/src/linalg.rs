use crate::error::{err, Result};

pub struct Csr {
    pub n: usize,
    pub indptr: Vec<usize>,
    pub indices: Vec<usize>,
    pub data: Vec<f64>,
}

pub fn csr_from_triplets(n: usize, mut trips: Vec<(usize, usize, f64)>) -> Csr {
    trips.retain(|(i, j, v)| *i < n && *j < n && v.abs() > 0.0);
    trips.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut indptr = vec![0usize; n + 1];
    let mut indices = Vec::with_capacity(trips.len());
    let mut data = Vec::with_capacity(trips.len());
    let mut k = 0;
    while k < trips.len() {
        let i = trips[k].0;
        let j = trips[k].1;
        let mut v = 0.0;
        while k < trips.len() && trips[k].0 == i && trips[k].1 == j {
            v += trips[k].2;
            k += 1;
        }
        indices.push(j);
        data.push(v);
        indptr[i + 1] += 1;
    }
    for i in 0..n {
        indptr[i + 1] += indptr[i];
    }
    Csr {
        n,
        indptr,
        indices,
        data,
    }
}

impl Csr {
    pub fn matvec(&self, x: &[f64], y: &mut [f64]) {
        y.fill(0.0);
        for i in 0..self.n {
            let mut s = 0.0;
            for k in self.indptr[i]..self.indptr[i + 1] {
                s += self.data[k] * x[self.indices[k]];
            }
            y[i] = s;
        }
    }

    pub fn diag(&self) -> Vec<f64> {
        let mut d = vec![0.0; self.n];
        for i in 0..self.n {
            for k in self.indptr[i]..self.indptr[i + 1] {
                if self.indices[k] == i {
                    d[i] = self.data[k];
                    break;
                }
            }
        }
        d
    }

    pub fn to_dense(&self) -> Vec<f64> {
        let n = self.n;
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            for k in self.indptr[i]..self.indptr[i + 1] {
                a[i * n + self.indices[k]] = self.data[k];
            }
        }
        a
    }
}

/// In-place Cholesky (lower) of dense SPD n×n, then solve A x = b.
pub fn chol_solve(a: &mut [f64], n: usize, b: &[f64]) -> Result<Vec<f64>> {
    for i in 0..n {
        for j in 0..=i {
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= a[i * n + k] * a[j * n + k];
            }
            if i == j {
                if s < 0.0 {
                    return err("Steifigkeitsmatrix ist indefinit.");
                }
                if s <= 1e-30 {
                    return err(
                        "Steifigkeitsmatrix ist singulär — Randbedingungen unzureichend (Starrkörperbewegung oder entartete Elemente).",
                    );
                }
                a[i * n + i] = s.sqrt();
            } else {
                a[i * n + j] = s / a[j * n + j];
            }
        }
    }
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= a[i * n + k] * y[k];
        }
        y[i] = s / a[i * n + i];
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for k in (i + 1)..n {
            s -= a[k * n + i] * x[k];
        }
        x[i] = s / a[i * n + i];
    }
    Ok(x)
}

/// Dense LU with partial pivoting. Symmetric indefinite tangents (Riks past
/// the limit point) are not SPD, so Cholesky refuses them.
pub fn lu_solve(a: &[f64], n: usize, b: &[f64]) -> Result<Vec<f64>> {
    if n == 0 || a.len() < n * n || b.len() < n {
        return err("LU: Dimension.");
    }
    let mut scale = 0.0_f64;
    for &v in a.iter().take(n * n) {
        scale = scale.max(v.abs());
    }
    if scale <= 0.0 {
        return err(
            "Steifigkeitsmatrix ist singulär — Randbedingungen unzureichend (Starrkörperbewegung oder entartete Elemente).",
        );
    }
    let mut lu = a[..n * n].to_vec();
    let mut piv: Vec<usize> = (0..n).collect();
    for k in 0..n {
        let mut pivrow = k;
        let mut maxv = lu[k * n + k].abs();
        for i in (k + 1)..n {
            let v = lu[i * n + k].abs();
            if v > maxv {
                maxv = v;
                pivrow = i;
            }
        }
        if maxv <= 1e-14 * scale {
            return err(
                "Steifigkeitsmatrix ist singulär — Randbedingungen unzureichend (Starrkörperbewegung oder entartete Elemente).",
            );
        }
        if pivrow != k {
            for j in 0..n {
                lu.swap(k * n + j, pivrow * n + j);
            }
            piv.swap(k, pivrow);
        }
        let akk = lu[k * n + k];
        for i in (k + 1)..n {
            lu[i * n + k] /= akk;
            let lik = lu[i * n + k];
            for j in (k + 1)..n {
                lu[i * n + j] -= lik * lu[k * n + j];
            }
        }
    }
    let mut y = vec![0.0; n];
    for i in 0..n {
        y[i] = b[piv[i]];
    }
    for i in 0..n {
        let mut s = y[i];
        for k in 0..i {
            s -= lu[i * n + k] * y[k];
        }
        y[i] = s;
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for k in (i + 1)..n {
            s -= lu[i * n + k] * x[k];
        }
        let diag = lu[i * n + i];
        if diag.abs() <= 1e-14 * scale {
            return err(
                "Steifigkeitsmatrix ist singulär — Randbedingungen unzureichend (Starrkörperbewegung oder entartete Elemente).",
            );
        }
        x[i] = s / diag;
    }
    Ok(x)
}

pub struct CgInfo {
    pub iters: usize,
    pub residual: f64,
}

pub fn pcg(a: &Csr, b: &[f64], tol: f64, max_iter: usize) -> Result<(Vec<f64>, CgInfo)> {
    let n = a.n;
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    let mut bnorm = 0.0;
    for &v in b {
        bnorm += v * v;
    }
    bnorm = bnorm.sqrt();
    if bnorm < 1e-30 {
        return Ok((
            x,
            CgInfo {
                iters: 0,
                residual: 0.0,
            },
        ));
    }
    let diag = a.diag();
    let mut z = vec![0.0; n];
    for i in 0..n {
        z[i] = r[i] / if diag[i].abs() > 1e-30 { diag[i] } else { 1.0 };
    }
    let mut p = z.clone();
    let mut rz: f64 = r.iter().zip(z.iter()).map(|(ri, zi)| ri * zi).sum();
    let mut ap = vec![0.0; n];
    let mut iters = 0;
    let mut residual = bnorm;
    for it in 0..max_iter {
        iters = it + 1;
        a.matvec(&p, &mut ap);
        let pap: f64 = p.iter().zip(ap.iter()).map(|(pi, ai)| pi * ai).sum();
        if pap.abs() < 1e-30 {
            return err("CG: unerwartetes Null-Pivot — Matrix nicht positiv definit.");
        }
        let alpha = rz / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        residual = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        if residual <= tol * bnorm {
            break;
        }
        for i in 0..n {
            z[i] = r[i] / if diag[i].abs() > 1e-30 { diag[i] } else { 1.0 };
        }
        let rz_new: f64 = r.iter().zip(z.iter()).map(|(ri, zi)| ri * zi).sum();
        let beta = rz_new / rz;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz = rz_new;
    }
    if residual > tol * bnorm * 10.0 {
        return err(format!(
            "CG konvergierte nicht (Residuum {residual:.3e} nach {iters} Iterationen)."
        ));
    }
    Ok((x, CgInfo { iters, residual }))
}

pub struct SparseResult {
    pub x: Vec<f64>,
    pub name: String,
    pub iters: usize,
    pub residual: f64,
}

fn residual_of(a: &Csr, x: &[f64], b: &[f64]) -> f64 {
    let mut ax = vec![0.0; a.n];
    a.matvec(x, &mut ax);
    let mut s = 0.0;
    for i in 0..a.n {
        let d = ax[i] - b[i];
        s += d * d;
    }
    s.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lu_solves_indefinite_2x2() {
        // diag(1, -2) * [3, 4] = [3, -8]
        let a = [1.0, 0.0, 0.0, -2.0];
        let x = lu_solve(&a, 2, &[3.0, -8.0]).unwrap();
        assert!((x[0] - 3.0).abs() < 1e-12);
        assert!((x[1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn cholesky_rejects_indefinite_and_lu_fallback_solves() {
        let trips = vec![(0, 0, 1.0), (1, 1, -2.0)];
        let rhs = [3.0, -8.0];
        let s = crate::backend::with_sparse_backend(crate::backend::SparseBackend::Cholesky, || {
            solve_kff(2, trips, &rhs)
        })
        .unwrap();
        assert_eq!(s.name, "dense LU");
        assert!((s.x[0] - 3.0).abs() < 1e-8, "{:?}", s.x);
        assert!((s.x[1] - 4.0).abs() < 1e-8, "{:?}", s.x);
    }
}

/// Factor and solve K_ff x = rhs.
///
/// Native: `--solver` / `AXIA_SOLVER` selects PARDISO, faer, rivrs-sparse,
/// dense Cholesky or PCG. Default `auto` is MKL → Panua → faer → rivrs.
/// WASM uses the in-crate Cholesky, and dense LU when the tangent is indefinite.
pub fn solve_kff(n: usize, trips: Vec<(usize, usize, f64)>, rhs: &[f64]) -> Result<SparseResult> {
    crate::backend::apply_env_solver();
    let want = crate::backend::sparse_backend();
    let csr = csr_from_triplets(n, trips);

    const DENSE_LIMIT: usize = 900;

    let dense = |csr: &Csr, forced: bool| -> Result<SparseResult> {
        if n > DENSE_LIMIT && forced {
            return err(format!(
                "dense Cholesky: n={n} > {DENSE_LIMIT}. Wähle --solver faer|rivrs|pcg."
            ));
        }
        if n > DENSE_LIMIT {
            return err("intern: dense Cholesky nur für kleine Systeme");
        }
        let mut a = csr.to_dense();
        match chol_solve(&mut a, n, rhs) {
            Ok(x) => {
                #[cfg(not(target_arch = "wasm32"))]
                eprintln!("axia: sparse solver: dense Cholesky");
                let residual = residual_of(csr, &x, rhs);
                Ok(SparseResult {
                    x,
                    name: "Cholesky".into(),
                    iters: 1,
                    residual,
                })
            }
            Err(e) => {
                if !e.to_string().contains("indefinit") {
                    return Err(e);
                }
                // Snap-through: K is symmetric indefinite. Faer uses LU here;
                // the in-crate path (WASM, `--solver cholesky`) must too.
                #[cfg(not(target_arch = "wasm32"))]
                eprintln!("axia: sparse solver: dense LU");
                let a = csr.to_dense();
                let x = lu_solve(&a, n, rhs)?;
                let residual = residual_of(csr, &x, rhs);
                Ok(SparseResult {
                    x,
                    name: "dense LU".into(),
                    iters: 1,
                    residual,
                })
            }
        }
    };
    let iterative = |csr: &Csr| -> Result<SparseResult> {
        #[cfg(not(target_arch = "wasm32"))]
        eprintln!("axia: sparse solver: PCG");
        let (x, info) = pcg(csr, rhs, 1e-8, (4 * n).max(200))?;
        Ok(SparseResult {
            x,
            name: "PCG".into(),
            iters: info.iters,
            residual: info.residual,
        })
    };

    if want == crate::backend::SparseBackend::Cholesky {
        return dense(&csr, true);
    }
    if want == crate::backend::SparseBackend::Pcg {
        return iterative(&csr);
    }

    #[cfg(not(target_arch = "wasm32"))]
    {
        match crate::sparse_native::solve_kff(&csr, rhs) {
            Ok(s) => {
                return Ok(SparseResult {
                    x: s.x,
                    name: s.name,
                    iters: s.iters,
                    residual: s.residual,
                });
            }
            Err(e) => {
                if want != crate::backend::SparseBackend::Auto {
                    return Err(e);
                }
                eprintln!("axia: sparse solver failed ({e}), falling back to in-crate solver");
            }
        }
    }

    if n <= DENSE_LIMIT {
        dense(&csr, false)
    } else {
        iterative(&csr)
    }
}
