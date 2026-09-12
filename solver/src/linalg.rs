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
