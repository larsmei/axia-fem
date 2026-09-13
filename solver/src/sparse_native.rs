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
    use std::sync::atomic::{AtomicBool, Ordering};
    static DONE: AtomicBool = AtomicBool::new(false);
    if !DONE.swap(true, Ordering::Relaxed) {
        eprintln!("axia: sparse solver: {name}");
    }
}

/// Names `pardiso-wrapper` 0.1.2 looks up (`libmkl_rt.dll` on Windows).
fn wrapper_lib_name() -> &'static str {
    if cfg!(target_os = "windows") {
        "libmkl_rt.dll"
    } else if cfg!(target_os = "macos") {
        "libmkl_rt.dylib"
    } else {
        "libmkl_rt.so"
    }
}

/// Intel oneAPI actually ships `mkl_rt.dll` / `mkl_rt.2.dll` (no `lib` prefix)
/// and `libmkl_rt.so[.2]` under `$MKLROOT/{bin,lib,…}`, not only `$MKLROOT/lib`.
fn mkl_candidate_names() -> &'static [&'static str] {
    // Search every known Intel filename so a Windows MKL tree is still
    // found when we unit-test the walker on Linux (and vice versa).
    &[
        "mkl_rt.dll",
        "mkl_rt.2.dll",
        "mkl_rt.1.dll",
        "libmkl_rt.dll",
        "libmkl_rt.so",
        "libmkl_rt.so.2",
        "libmkl_rt.so.1",
        "libmkl_rt.so.2.0",
        "libmkl_rt.dylib",
        "libmkl_rt.2.dylib",
        "libmkl_rt.1.dylib",
    ]
}

fn mkl_subdirs() -> &'static [&'static str] {
    &[
        "bin",
        "bin/intel64",
        "redist/intel64",
        "redist/intel64/msmpi",
        "lib",
        "lib/intel64",
        "lib/intel64_win",
        "lib/intel64/lib",
        "",
    ]
}

fn path_prepend(dir: &std::path::Path) {
    let key = if cfg!(target_os = "windows") {
        "PATH"
    } else if cfg!(target_os = "macos") {
        "DYLD_LIBRARY_PATH"
    } else {
        "LD_LIBRARY_PATH"
    };
    let dir_s = dir.to_string_lossy().into_owned();
    let joined = match std::env::var_os(key) {
        Some(old) => {
            let mut v = vec![std::ffi::OsString::from(&dir_s)];
            v.push(old);
            std::env::join_paths(v).unwrap_or_else(|_| std::ffi::OsString::from(&dir_s))
        }
        None => std::ffi::OsString::from(&dir_s),
    };
    unsafe {
        std::env::set_var(key, joined);
    }
}

fn compiler_redist_dirs(mkl_bin: &std::path::Path) -> Vec<std::path::PathBuf> {
    let mut out = Vec::new();
    let mut push = |p: std::path::PathBuf| {
        if p.is_dir() && !out.contains(&p) {
            out.push(p);
        }
    };
    if let Some(oneapi) = std::env::var_os("ONEAPI_ROOT") {
        let root = std::path::PathBuf::from(oneapi);
        push(root.join("compiler/latest/bin"));
        push(root.join("compiler/latest/windows/redist/intel64_win"));
        push(root.join("compiler/latest/lib"));
    }
    // $MKLROOT/bin → ../../compiler/latest/bin
    if let Some(mkl_root) = mkl_bin.parent() {
        push(mkl_root.join("../compiler/latest/bin"));
        push(mkl_root.join("../../compiler/latest/bin"));
        push(mkl_root.join("../../compiler/latest/windows/redist/intel64_win"));
        push(mkl_root.join("../compiler/latest/lib"));
    }
    out
}

/// Locate a real MKL runtime library on disk (does not mutate env).
pub(crate) fn find_mkl_runtime() -> Option<std::path::PathBuf> {
    find_mkl_runtime_in(mkl_search_roots())
}

fn mkl_search_roots() -> Vec<std::path::PathBuf> {
    let mut roots = Vec::new();
    let mut push = |p: std::path::PathBuf| {
        if !p.as_os_str().is_empty() && !roots.contains(&p) {
            roots.push(p);
        }
    };
    if let Some(p) = std::env::var_os("MKL_PARDISO_PATH") {
        let pb = std::path::PathBuf::from(p);
        if pb.is_file() {
            if let Some(d) = pb.parent() {
                push(d.to_path_buf());
            }
        } else {
            push(pb);
        }
    }
    if let Some(p) = std::env::var_os("MKLROOT") {
        push(std::path::PathBuf::from(p));
    }
    if let Some(p) = std::env::var_os("ONEAPI_ROOT") {
        let root = std::path::PathBuf::from(p);
        push(root.join("mkl/latest"));
        push(root.join("mkl"));
    }
    push(std::path::PathBuf::from("/opt/intel/oneapi/mkl/latest"));
    push(std::path::PathBuf::from("/opt/intel/mkl"));
    roots
}

fn find_mkl_runtime_in(roots: Vec<std::path::PathBuf>) -> Option<std::path::PathBuf> {
    let names = mkl_candidate_names();
    for root in &roots {
        for sub in mkl_subdirs() {
            let dir = if sub.is_empty() {
                root.clone()
            } else {
                root.join(sub)
            };
            for name in names {
                let p = dir.join(name);
                if p.is_file() {
                    return Some(p);
                }
            }
        }
    }
    // PATH / LD_LIBRARY_PATH entries (setvars.bat puts MKL bin here)
    let path_key = if cfg!(target_os = "windows") {
        "PATH"
    } else if cfg!(target_os = "macos") {
        "DYLD_LIBRARY_PATH"
    } else {
        "LD_LIBRARY_PATH"
    };
    if let Some(paths) = std::env::var_os(path_key) {
        for dir in std::env::split_paths(&paths) {
            for name in names {
                let p = dir.join(name);
                if p.is_file() {
                    return Some(p);
                }
            }
        }
    }
    None
}

/// `pardiso-wrapper` 0.1.2 uses `lazy_static` and looks for `libmkl_rt.dll` in
/// `$MKLROOT/lib` only. Intel oneAPI ships `mkl_rt.dll` under `$MKLROOT/bin`.
/// Must run **before** the first `is_available()` call.
fn prepare_mkl_env() {
    use std::sync::Once;
    static ONCE: Once = Once::new();
    ONCE.call_once(prepare_mkl_env_inner);
}

fn prepare_mkl_env_inner() {
    let expected = wrapper_lib_name();
    let Some(found) = find_mkl_runtime() else {
        if std::env::var_os("MKLROOT").is_some() {
            eprintln!(
                "axia: MKLROOT is set but no MKL runtime was found \
                 (looked for mkl_rt / libmkl_rt under $MKLROOT/bin, lib, redist). \
                 Falling back to rivrs-sparse. Hint: run oneAPI setvars, or put \
                 the folder that contains mkl_rt.dll on PATH."
            );
        }
        return;
    };
    let Some(dir) = found.parent().map(|p| p.to_path_buf()) else {
        return;
    };

    path_prepend(&dir);
    for extra in compiler_redist_dirs(&dir) {
        path_prepend(&extra);
    }

    let name = found.file_name().and_then(|s| s.to_str()).unwrap_or("");
    if name == expected {
        if std::env::var_os("MKL_PARDISO_PATH").is_none() {
            unsafe {
                std::env::set_var("MKL_PARDISO_PATH", &dir);
            }
        }
        return;
    }

    // Wrapper only searches for `libmkl_rt.dll` — copy/hardlink under that name.
    let shim_dir = std::env::temp_dir().join("axia-mkl-shim");
    if let Err(e) = std::fs::create_dir_all(&shim_dir) {
        eprintln!("axia: cannot create MKL shim dir {}: {e}", shim_dir.display());
        return;
    }
    let shim = shim_dir.join(expected);
    if !shim.exists() {
        let linked = std::fs::hard_link(&found, &shim).is_ok();
        if !linked {
            if let Err(e) = std::fs::copy(&found, &shim) {
                eprintln!(
                    "axia: cannot create MKL shim {} → {}: {e}",
                    found.display(),
                    shim.display()
                );
                return;
            }
        }
    }
    path_prepend(&shim_dir);
    unsafe {
        std::env::set_var("MKL_PARDISO_PATH", &shim_dir);
    }
    eprintln!(
        "axia: MKL runtime {} (wrapper expects {expected}; shim {})",
        found.display(),
        shim.display()
    );
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
    prepare_mkl_env();
    #[cfg(target_arch = "x86_64")]
    {
        if pardiso_wrapper::MKLPardisoSolver::is_available() {
            let name = "PARDISO (Intel MKL)";
            announce(name);
            match run_pardiso::<pardiso_wrapper::MKLPardisoSolver>(csr, rhs) {
                Ok(x) => return Some((x, name.to_string())),
                Err(e) => eprintln!("axia: {name} failed ({e}), falling back"),
            }
        } else if std::env::var_os("MKLROOT").is_some() {
            use std::sync::atomic::{AtomicBool, Ordering};
            static HINT: AtomicBool = AtomicBool::new(false);
            if !HINT.swap(true, Ordering::Relaxed) {
                eprintln!(
                    "axia: MKLROOT is set but PARDISO (Intel MKL) did not load. \
                     Need mkl_rt.dll (shimmed as libmkl_rt.dll) and OpenMP \
                     (libiomp5md.dll, typically oneAPI compiler/latest/bin) on PATH. \
                     Falling back to rivrs-sparse."
                );
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_intel_windows_dll_under_bin() {
        let tmp = std::env::temp_dir().join(format!(
            "axia-mkl-find-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let bin = tmp.join("bin");
        std::fs::create_dir_all(&bin).unwrap();
        let dll = bin.join("mkl_rt.dll");
        std::fs::write(&dll, b"fake-mkl").unwrap();
        let found = find_mkl_runtime_in(vec![tmp.clone()]).expect("should find mkl_rt.dll");
        assert_eq!(found, dll);
        // libmkl_rt.dll (wrapper name) is NOT required to exist
        assert!(!bin.join("libmkl_rt.dll").exists());
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn finds_so_dot_two_under_lib() {
        let tmp = std::env::temp_dir().join(format!(
            "axia-mkl-so-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let lib = tmp.join("lib");
        std::fs::create_dir_all(&lib).unwrap();
        let so = lib.join("libmkl_rt.so.2");
        std::fs::write(&so, b"fake-mkl").unwrap();
        let found = find_mkl_runtime_in(vec![tmp.clone()]).expect("should find libmkl_rt.so.2");
        assert_eq!(found, so);
        let _ = std::fs::remove_dir_all(&tmp);
    }
}
