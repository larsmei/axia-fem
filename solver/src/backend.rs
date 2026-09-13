//! Sparse solver selection (CLI `--solver` / env `AXIA_SOLVER`).

use std::cell::Cell;
use std::sync::atomic::{AtomicU8, Ordering};

use crate::error::{err, Result};

/// Native sparse / dense backends that can factor \(K_{ff}x=b\).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum SparseBackend {
    /// MKL → Panua → faer → rivrs → dense/PCG.
    Auto = 0,
    /// Intel MKL PARDISO only (needs `mkl_rt`).
    Mkl = 1,
    /// Panua PARDISO only.
    Panua = 2,
    /// Either PARDISO flavour, no pure-Rust fallback.
    Pardiso = 3,
    /// faer supernodal \(LL^\top\) (then \(LU\)). Pure Rust, no extra libs.
    Faer = 4,
    /// rivrs-sparse multifrontal \(LDL^\top\) (APTP). Pure Rust (+ optional METIS).
    Rivrs = 5,
    /// In-crate dense Cholesky (small \(n\)).
    Cholesky = 6,
    /// In-crate PCG.
    Pcg = 7,
}

impl SparseBackend {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Mkl => "mkl",
            Self::Panua => "panua",
            Self::Pardiso => "pardiso",
            Self::Faer => "faer",
            Self::Rivrs => "rivrs",
            Self::Cholesky => "cholesky",
            Self::Pcg => "pcg",
        }
    }

    pub fn is_dense_path(self) -> bool {
        matches!(self, Self::Cholesky | Self::Pcg)
    }
}

static PROCESS: AtomicU8 = AtomicU8::new(SparseBackend::Auto as u8);

thread_local! {
    static TLS: Cell<Option<SparseBackend>> = const { Cell::new(None) };
}

/// Process-wide default (CLI `--solver`, env `AXIA_SOLVER`).
pub fn set_sparse_backend(b: SparseBackend) {
    PROCESS.store(b as u8, Ordering::Relaxed);
}

pub fn sparse_backend() -> SparseBackend {
    if let Some(b) = TLS.with(|t| t.get()) {
        return b;
    }
    from_u8(PROCESS.load(Ordering::Relaxed))
}

/// Apply `AXIA_SOLVER` once if the process default is still `auto` and the
/// environment names a backend.
pub fn apply_env_solver() {
    if PROCESS.load(Ordering::Relaxed) != SparseBackend::Auto as u8 {
        return;
    }
    let Ok(raw) = std::env::var("AXIA_SOLVER") else {
        return;
    };
    if let Ok(b) = parse_sparse_backend(&raw) {
        PROCESS.store(b as u8, Ordering::Relaxed);
    }
}

/// Run `f` with a thread-local backend (tests). Restores the previous value.
pub fn with_sparse_backend<R>(b: SparseBackend, f: impl FnOnce() -> R) -> R {
    TLS.with(|t| {
        let prev = t.replace(Some(b));
        let r = f();
        t.set(prev);
        r
    })
}

pub fn parse_sparse_backend(s: &str) -> Result<SparseBackend> {
    let t = s.trim().to_ascii_lowercase();
    let t = t.strip_prefix("--solver=").unwrap_or(&t);
    Ok(match t {
        "auto" | "default" => SparseBackend::Auto,
        "mkl" | "intel" | "pardiso-mkl" | "intel-mkl" => SparseBackend::Mkl,
        "panua" | "pardiso-panua" => SparseBackend::Panua,
        "pardiso" => SparseBackend::Pardiso,
        "faer" | "supernodal" | "llt" => SparseBackend::Faer,
        "rivrs" | "rivrs-sparse" | "ldlt" => SparseBackend::Rivrs,
        "cholesky" | "dense" | "chol" => SparseBackend::Cholesky,
        "pcg" | "cg" => SparseBackend::Pcg,
        _ => {
            return err(format!(
                "unbekannter Sparse-Solver '{s}'. \
                 Erlaubt: auto, mkl, panua, pardiso, faer, rivrs, cholesky, pcg"
            ));
        }
    })
}

fn from_u8(v: u8) -> SparseBackend {
    match v {
        1 => SparseBackend::Mkl,
        2 => SparseBackend::Panua,
        3 => SparseBackend::Pardiso,
        4 => SparseBackend::Faer,
        5 => SparseBackend::Rivrs,
        6 => SparseBackend::Cholesky,
        7 => SparseBackend::Pcg,
        _ => SparseBackend::Auto,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_aliases() {
        assert_eq!(parse_sparse_backend("FAER").unwrap(), SparseBackend::Faer);
        assert_eq!(parse_sparse_backend("rivrs-sparse").unwrap(), SparseBackend::Rivrs);
        assert_eq!(parse_sparse_backend("intel").unwrap(), SparseBackend::Mkl);
        assert!(parse_sparse_backend("umfpack").is_err());
    }
}
