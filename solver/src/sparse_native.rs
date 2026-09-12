//! Native sparse backends: PARDISO (MKL / Panua) with rivrs-sparse as fallback.

use pardiso_wrapper::{MessageLevel, MatrixType, PardisoInterface, Phase};

use crate::error::{err, FemError, Result};
use crate::linalg::Csr;

pub struct SparseSolve {
    pub x: Vec<f64>,
    pub name: String,
    pub iters: usize,
    pub residual: f64,
}

fn residual(csr: &Csr, x: &[f64], b: &[f64]) -> f64 {
    let mut ax = vec![0.0; csr.n];
    csr.matvec(x, &mut ax);
    let mut s = 0.0;
    for i in 0..csr.n {
        let d = ax[i] - b[i];
        s += d * d;
    }
    s.sqrt()
}

fn announce(name: &str) {
    eprintln!("axia: sparse solver: {name}");
}

/// 1-based CSR of the upper triangle (incl. diagonal), columns sorted per row.
fn csr_upper_1based(csr: &Csr) -> Result<(Vec<f64>, Vec<i32>, Vec<i32>)> {
    if csr.n > i32::MAX as usize {
        return err("System zu groß für PARDISO (n > i32::MAX).");
    }
    let n = csr.n;
    let mut ia = vec![0i32; n + 1];
    let mut ja = Vec::new();
    let mut a = Vec::new();
    for i in 0..n {
        ia[i] = (ja.len() as i32) + 1;
        let mut cols: Vec<(usize, f64)> = Vec::new();
        for k in csr.indptr[i]..csr.indptr[i + 1] {
            let j = csr.indices[k];
            if j >= i {
                cols.push((j, csr.data[k]));
            }
        }
        cols.sort_by_key(|c| c.0);
        // merge duplicates
        let mut m: Vec<(usize, f64)> = Vec::new();
        for (j, v) in cols {
            if let Some(last) = m.last_mut() {
                if last.0 == j {
                    last.1 += v;
                    continue;
                }
            }
            m.push((j, v));
        }
        if m.iter().all(|(j, _)| *j != i) {
            m.push((i, 0.0));
            m.sort_by_key(|c| c.0);
        }
        for (j, v) in m {
            ja.push(j as i32 + 1);
            a.push(v);
        }
    }
    ia[n] = ja.len() as i32 + 1;
    Ok((a, ia, ja))
}

fn run_pardiso<S: PardisoInterface>(csr: &Csr, rhs: &[f64]) -> Result<Vec<f64>> {
    let n = csr.n as i32;
    let (a, ia, ja) = csr_upper_1based(csr)?;
    let try_type = |mtype: MatrixType| -> Result<Vec<f64>> {
        let mut b = rhs.to_vec();
        let mut x = vec![0.0; csr.n];
        let mut ps = S::new().map_err(|e| FemError(e.to_string()))?;
        ps.set_matrix_type(mtype);
        ps.pardisoinit().map_err(|e| FemError(e.to_string()))?;
        ps.set_message_level(MessageLevel::Off);
        ps.set_phase(Phase::AnalysisNumFactSolveRefine);
        ps.pardiso(&a, &ia, &ja, &mut b, &mut x, n, 1)
            .map_err(|e| FemError(e.to_string()))?;
        Ok(x)
    };
    match try_type(MatrixType::RealSymmetricPositiveDefinite) {
        Ok(x) => Ok(x),
        Err(_) => try_type(MatrixType::RealSymmetricIndefinite),
    }
}

fn try_pardiso(csr: &Csr, rhs: &[f64]) -> Option<(Vec<f64>, String)> {
    #[cfg(target_arch = "x86_64")]
    {
        if pardiso_wrapper::MKLPardisoSolver::is_available() {
            let name = "PARDISO (Intel MKL)";
            announce(name);
            match run_pardiso::<pardiso_wrapper::MKLPardisoSolver>(csr, rhs) {
                Ok(x) => return Some((x, name.to_string())),
                Err(e) => eprintln!("axia: {name} failed ({e}), falling back"),
            }
        }
    }
    if pardiso_wrapper::PanuaPardisoSolver::is_available() {
        let name = "PARDISO (Panua)";
        announce(name);
        match run_pardiso::<pardiso_wrapper::PanuaPardisoSolver>(csr, rhs) {
            Ok(x) => return Some((x, name.to_string())),
            Err(e) => eprintln!("axia: {name} failed ({e}), falling back"),
        }
    }
    None
}

fn csr_to_faer(csr: &Csr) -> Result<faer::sparse::SparseColMat<usize, f64>> {
    use faer::sparse::{SparseColMat, Triplet};
    let mut trips = Vec::with_capacity(csr.data.len());
    for i in 0..csr.n {
        for k in csr.indptr[i]..csr.indptr[i + 1] {
            trips.push(Triplet::new(i, csr.indices[k], csr.data[k]));
        }
    }
    SparseColMat::try_new_from_triplets(csr.n, csr.n, &trips)
        .map_err(|e| FemError(format!("rivrs-sparse: Matrixaufbau fehlgeschlagen ({e})")))
}

fn try_rivrs(csr: &Csr, rhs: &[f64]) -> Result<Vec<f64>> {
    use faer::{Col, Par};
    use rivrs_sparse::symmetric::{OrderingStrategy, SolverOptions, SparseLDLT};

    let mat = csr_to_faer(csr)?;
    let b = Col::from_fn(csr.n, |i| rhs[i]);
    let orderings = [OrderingStrategy::Amd, OrderingStrategy::Metis];
    let mut last = String::new();
    for ordering in &orderings {
        let tag = match ordering {
            OrderingStrategy::Amd => "AMD",
            OrderingStrategy::Metis => "METIS",
            _ => "custom",
        };
        let mut opts = SolverOptions::default();
        opts.ordering = ordering.clone();
        opts.par = Par::Seq;
        match SparseLDLT::solve_full(&mat, &b, &opts) {
            Ok(x) => {
                return Ok((0..csr.n).map(|i| x[i]).collect());
            }
            Err(e) => last = format!("{tag}: {e}"),
        }
    }
    err(format!("rivrs-sparse failed ({last})"))
}

pub fn solve_kff(csr: &Csr, rhs: &[f64]) -> Result<SparseSolve> {
    if let Some((x, name)) = try_pardiso(csr, rhs) {
        let residual = residual(csr, &x, rhs);
        return Ok(SparseSolve {
            x,
            name,
            iters: 1,
            residual,
        });
    }

    let name = "rivrs-sparse (LDLT)";
    announce(name);
    match try_rivrs(csr, rhs) {
        Ok(x) => {
            let residual = residual(csr, &x, rhs);
            Ok(SparseSolve {
                x,
                name: name.to_string(),
                iters: 1,
                residual,
            })
        }
        Err(e) => Err(e),
    }
}
