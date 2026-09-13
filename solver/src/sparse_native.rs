//! Native sparse backends: PARDISO (MKL / Panua) with rivrs-sparse as fallback.

use pardiso_wrapper::{MessageLevel, MatrixType, PardisoInterface, Phase};

use crate::backend::SparseBackend;
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
    use std::sync::Mutex;
    static LAST: Mutex<Option<String>> = Mutex::new(None);
    if let Ok(mut g) = LAST.lock() {
        if g.as_deref() != Some(name) {
            eprintln!("axia: sparse solver: {name}");
            *g = Some(name.to_string());
        }
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
        "", // DLLs dropped next to axia.exe / in cwd
        "bin",
        "bin/intel64",
        "redist/intel64",
        "redist/intel64/msmpi",
        "lib",
        "lib/intel64",
        "lib/intel64_win",
        "lib/intel64/lib",
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
    if let Some(mkl_root) = mkl_bin.parent() {
        push(mkl_root.join("../compiler/latest/bin"));
        push(mkl_root.join("../../compiler/latest/bin"));
        push(mkl_root.join("../../compiler/latest/windows/redist/intel64_win"));
        push(mkl_root.join("../compiler/latest/lib"));
    }
    out
}

fn exe_dir() -> Option<std::path::PathBuf> {
    std::env::current_exe()
        .ok()
        .and_then(|p| p.parent().map(|d| d.to_path_buf()))
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
    // Portable Windows layout: mkl_rt.dll next to axia.exe, even if cwd differs.
    if let Some(d) = exe_dir() {
        push(d);
    }
    if let Ok(d) = std::env::current_dir() {
        push(d);
    }
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

fn is_mkl_companion_name(name: &str) -> bool {
    let n = name.to_ascii_lowercase();
    let ext_ok = n.ends_with(".dll")
        || n.ends_with(".so")
        || n.contains(".so.")
        || n.ends_with(".dylib");
    if !ext_ok {
        return false;
    }
    n.starts_with("mkl_")
        || n.starts_with("libmkl_")
        || n.contains("iomp")
        || n.contains("libomp")
        || n.starts_with("libiomp")
}

fn is_lib_filename(name: &str) -> bool {
    let n = name.to_ascii_lowercase();
    n.ends_with(".dll") || n.ends_with(".so") || n.contains(".so.") || n.ends_with(".dylib")
}

fn is_mkl_core_name(name: &str) -> bool {
    is_lib_filename(name) && name.to_ascii_lowercase().contains("mkl_core")
}

fn is_mkl_thread_layer_name(name: &str) -> bool {
    if !is_lib_filename(name) {
        return false;
    }
    let n = name.to_ascii_lowercase();
    n.contains("mkl_sequential")
        || n.contains("mkl_intel_thread")
        || n.contains("mkl_tbb_thread")
        || n.contains("mkl_gnu_thread")
}

fn is_mkl_cpu_kernel_name(name: &str) -> bool {
    if !is_lib_filename(name) {
        return false;
    }
    let n = name.to_ascii_lowercase();
    let s = n.trim_start_matches("lib");
    s.starts_with("mkl_avx")
        || s.starts_with("mkl_def")
        || s.starts_with("mkl_mc")
        || s.starts_with("mkl_sse")
        || s.starts_with("mkl_p4")
}

fn is_iomp_name(name: &str) -> bool {
    if !is_lib_filename(name) {
        return false;
    }
    let n = name.to_ascii_lowercase();
    n.contains("iomp5") || n.starts_with("libiomp")
}

fn has_file_matching(dirs: &[std::path::PathBuf], pred: impl Fn(&str) -> bool) -> bool {
    for d in dirs {
        let Ok(rd) = std::fs::read_dir(d) else {
            continue;
        };
        for e in rd.flatten() {
            let name = e.file_name();
            let Some(s) = name.to_str() else { continue };
            if pred(s) && e.path().is_file() {
                return true;
            }
        }
    }
    false
}

fn missing_mkl_parts(dirs: &[std::path::PathBuf]) -> Vec<&'static str> {
    let mut missing = Vec::new();
    if !has_file_matching(dirs, is_mkl_core_name) {
        missing.push("mkl_core");
    }
    if !has_file_matching(dirs, is_mkl_thread_layer_name) {
        missing.push("mkl_sequential/mkl_intel_thread");
    }
    if !has_file_matching(dirs, is_mkl_cpu_kernel_name) {
        missing.push("mkl_avx2/mkl_def (CPU kernel)");
    }
    missing
}

fn mkl_source_dirs(rt_dir: &std::path::Path) -> Vec<std::path::PathBuf> {
    let mut dirs = Vec::new();
    let mut push = |p: std::path::PathBuf| {
        if p.is_dir() && !dirs.contains(&p) {
            dirs.push(p);
        }
    };
    push(rt_dir.to_path_buf());
    if let Some(d) = exe_dir() {
        push(d);
    }
    if let Ok(d) = std::env::current_dir() {
        push(d);
    }
    for extra in compiler_redist_dirs(rt_dir) {
        push(extra);
    }
    if let Some(parent) = rt_dir.parent() {
        push(parent.to_path_buf());
        push(parent.join("bin"));
        push(parent.join("lib"));
        push(parent.join("redist/intel64"));
        if let Some(gp) = parent.parent() {
            push(gp.join("bin"));
            push(gp.join("lib"));
            push(gp.join("compiler/latest/bin"));
            push(gp.join("compiler/latest/lib"));
        }
    }
    dirs
}

fn infer_mklroot(rt: &std::path::Path) -> Option<std::path::PathBuf> {
    let mut dir = rt.parent()?.to_path_buf();
    for _ in 0..3 {
        let name = dir
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        if matches!(
            name.as_str(),
            "bin" | "lib" | "intel64" | "intel64_win" | "redist"
        ) {
            dir = dir.parent()?.to_path_buf();
            continue;
        }
        break;
    }
    Some(dir)
}

fn choose_mkl_threading(dirs: &[std::path::PathBuf]) -> &'static str {
    if let Ok(v) = std::env::var("MKL_THREADING_LAYER") {
        if !v.trim().is_empty() {
            return "keep";
        }
    }
    let intel = has_file_matching(dirs, |n| n.to_ascii_lowercase().contains("mkl_intel_thread"));
    let seq = has_file_matching(dirs, |n| n.to_ascii_lowercase().contains("mkl_sequential"));
    let iomp = has_file_matching(dirs, is_iomp_name);
    if intel && iomp {
        "INTEL"
    } else if seq {
        "SEQUENTIAL"
    } else if intel {
        "INTEL"
    } else {
        "SEQUENTIAL"
    }
}

fn populate_bundle(dest: &std::path::Path, sources: &[std::path::PathBuf]) {
    let _ = std::fs::create_dir_all(dest);
    for src in sources {
        if src == dest {
            continue;
        }
        copy_mkl_companions(src, dest);
    }
}

fn copy_mkl_companions(from: &std::path::Path, to: &std::path::Path) {
    let Ok(rd) = std::fs::read_dir(from) else {
        return;
    };
    for e in rd.flatten() {
        let name = e.file_name();
        let Some(s) = name.to_str() else { continue };
        if !is_mkl_companion_name(s) {
            continue;
        }
        let dest = to.join(&name);
        if dest.exists() {
            continue;
        }
        let _ = std::fs::copy(e.path(), dest);
    }
}

/// `mkl_rt` loads `mkl_core` / threading layers from **its own directory**.
/// The wrapper name (`libmkl_rt.dll`) must therefore live next to those files,
/// not in a temp folder that only contains the shim.
fn ensure_wrapper_named_library(found: &std::path::Path) -> std::io::Result<std::path::PathBuf> {
    let expected = wrapper_lib_name();
    let name = found.file_name().and_then(|s| s.to_str()).unwrap_or("");
    if name.eq_ignore_ascii_case(expected) {
        return Ok(found.to_path_buf());
    }
    let dir = found.parent().ok_or_else(|| {
        std::io::Error::new(std::io::ErrorKind::NotFound, "MKL path has no parent")
    })?;
    let sibling = dir.join(expected);
    if sibling.exists() {
        return Ok(sibling);
    }
    if std::fs::hard_link(found, &sibling).is_ok() || std::fs::copy(found, &sibling).is_ok() {
        return Ok(sibling);
    }
    // Directory not writable (typical: Program Files). Copy the runtime plus
    // every MKL/OpenMP DLL from that folder into a user-writable shim dir so
    // mkl_rt can still find mkl_core.2.dll next to itself.
    let shim_dir = std::env::temp_dir().join("axia-mkl-shim");
    std::fs::create_dir_all(&shim_dir)?;
    let dest = shim_dir.join(expected);
    if !dest.exists() {
        if std::fs::hard_link(found, &dest).is_err() {
            std::fs::copy(found, &dest)?;
        }
    }
    copy_mkl_companions(dir, &shim_dir);
    Ok(dest)
}

#[cfg(windows)]
fn set_dll_directory(dir: &std::path::Path) {
    use std::os::windows::ffi::OsStrExt;
    let wide: Vec<u16> = dir.as_os_str().encode_wide().chain(Some(0)).collect();
    extern "system" {
        fn SetDllDirectoryW(lp_path_name: *const u16) -> i32;
    }
    unsafe {
        SetDllDirectoryW(wide.as_ptr());
    }
}

#[cfg(not(windows))]
fn set_dll_directory(_dir: &std::path::Path) {}

#[cfg(windows)]
fn probe_load(path: &std::path::Path) -> String {
    use std::os::windows::ffi::OsStrExt;
    let wide: Vec<u16> = path.as_os_str().encode_wide().chain(Some(0)).collect();
    extern "system" {
        fn LoadLibraryW(name: *const u16) -> isize;
        fn GetLastError() -> u32;
        fn FreeLibrary(h: isize) -> i32;
    }
    unsafe {
        let h = LoadLibraryW(wide.as_ptr());
        if h == 0 {
            let err = GetLastError();
            let hint = match err {
                126 => " (ERROR_MOD_NOT_FOUND: missing mkl_core.2.dll / mkl_intel_thread.2.dll / libiomp5md.dll next to mkl_rt.dll)",
                193 => " (not a valid Win32 image — 32/64-bit mismatch?)",
                _ => "",
            };
            format!("LoadLibrary({}) failed, Win32 {err}{hint}", path.display())
        } else {
            FreeLibrary(h);
            format!("LoadLibrary({}) ok", path.display())
        }
    }
}

#[cfg(not(windows))]
fn probe_load(path: &std::path::Path) -> String {
    format!("found {}", path.display())
}

fn mkl_dir_inventory(dir: &std::path::Path) -> String {
    let markers = [
        "mkl_rt.dll",
        "mkl_rt.2.dll",
        "libmkl_rt.dll",
        "mkl_core.2.dll",
        "mkl_core.dll",
        "mkl_intel_thread.2.dll",
        "mkl_sequential.2.dll",
        "mkl_avx2.2.dll",
        "mkl_def.2.dll",
        "libiomp5md.dll",
        "libmkl_rt.so",
        "libmkl_core.so.2",
        "libmkl_avx2.so.2",
        "libmkl_def.so.2",
        "libiomp5.so",
    ];
    let mut have = Vec::new();
    let mut missing = Vec::new();
    for n in markers {
        if dir.join(n).is_file() {
            have.push(n);
        } else if n.ends_with(".dll") {
            missing.push(n);
        }
    }
    format!(
        "{}: have [{}]{}",
        dir.display(),
        have.join(", "),
        if missing.is_empty() {
            String::new()
        } else {
            format!("; not found [{}]", missing.join(", "))
        }
    )
}

use std::sync::Mutex;
static MKL_DIAG: Mutex<Option<String>> = Mutex::new(None);
static MKL_BLOCKED: Mutex<Option<String>> = Mutex::new(None);

fn set_mkl_diag(s: impl Into<String>) {
    if let Ok(mut g) = MKL_DIAG.lock() {
        *g = Some(s.into());
    }
}

fn peek_mkl_diag() -> Option<String> {
    MKL_DIAG.lock().ok().and_then(|g| g.clone())
}

fn mkl_fail_detail() -> String {
    mkl_blocked_msg()
        .or_else(peek_mkl_diag)
        .unwrap_or_default()
}

fn set_mkl_blocked(s: impl Into<String>) {
    let s = s.into();
    set_mkl_diag(s.clone());
    if let Ok(mut g) = MKL_BLOCKED.lock() {
        *g = Some(s);
    }
}

fn mkl_blocked_msg() -> Option<String> {
    MKL_BLOCKED.lock().ok().and_then(|g| g.clone())
}

fn apply_mkl_env_defaults(dirs: &[std::path::PathBuf]) {
    unsafe {
        if std::env::var_os("KMP_DUPLICATE_LIB_OK").is_none() {
            // Mecway / other hosts often already loaded libomp. Without this,
            // Intel OpenMP prints Error #15 and abort()s — no Rust error.
            std::env::set_var("KMP_DUPLICATE_LIB_OK", "TRUE");
        }
        if std::env::var_os("MKL_INTERFACE_LAYER").is_none() {
            std::env::set_var("MKL_INTERFACE_LAYER", "LP64");
        }
        let layer = choose_mkl_threading(dirs);
        if layer != "keep" && std::env::var_os("MKL_THREADING_LAYER").is_none() {
            std::env::set_var("MKL_THREADING_LAYER", layer);
        }
    }
}

/// `pardiso-wrapper` 0.1.2 uses `lazy_static` and looks for `libmkl_rt.dll` in
/// `$MKLROOT/lib` only. Intel oneAPI ships `mkl_rt.dll` under `$MKLROOT/bin`
/// (and users drop it next to `axia.exe`). Must run **before** `is_available()`.
fn prepare_mkl_env() {
    use std::sync::Once;
    static ONCE: Once = Once::new();
    ONCE.call_once(prepare_mkl_env_inner);
}

fn prepare_mkl_env_inner() {
    let expected = wrapper_lib_name();

    // Always put the exe folder and cwd on the loader path first so a
    // portable copy of libmkl_rt.dll next to axia.exe is visible to `which`.
    if let Some(d) = exe_dir() {
        path_prepend(&d);
        set_dll_directory(&d);
    }
    if let Ok(d) = std::env::current_dir() {
        path_prepend(&d);
    }

    let Some(found) = find_mkl_runtime() else {
        let where_ = exe_dir()
            .map(|d| d.display().to_string())
            .unwrap_or_else(|| "<exe dir unknown>".into());
        set_mkl_diag(format!(
            "no mkl_rt / libmkl_rt next to axia.exe ({where_}), in cwd, $MKLROOT or PATH. \
             Wrapper looks for {expected}. Intel ships mkl_rt.dll. \
             Copy the MKL redist (mkl_rt.dll, mkl_core.2.dll, mkl_intel_thread.2.dll, \
             mkl_avx2.2.dll or mkl_def.2.dll, libiomp5md.dll) next to axia.exe, or set MKLROOT."
        ));
        return;
    };
    let Some(dir) = found.parent().map(|p| p.to_path_buf()) else {
        return;
    };

    path_prepend(&dir);
    set_dll_directory(&dir);
    for extra in compiler_redist_dirs(&dir) {
        path_prepend(&extra);
    }

    if std::env::var_os("MKLROOT").is_none() {
        if let Some(root) = infer_mklroot(&found) {
            unsafe {
                std::env::set_var("MKLROOT", &root);
            }
        }
    }

    let sources = mkl_source_dirs(&dir);
    apply_mkl_env_defaults(&sources);

    match ensure_wrapper_named_library(&found) {
        Ok(shim) => {
            if let Some(shim_dir) = shim.parent() {
                populate_bundle(shim_dir, &sources);
                path_prepend(shim_dir);
                set_dll_directory(shim_dir);
                unsafe {
                    std::env::set_var("MKL_PARDISO_PATH", shim_dir);
                }
            }
            let bundle = shim
                .parent()
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| dir.clone());
            let mut scan = sources.clone();
            if !scan.contains(&bundle) {
                scan.insert(0, bundle.clone());
            }
            let missing = missing_mkl_parts(&scan);
            let inv = mkl_dir_inventory(&dir);
            if !missing.is_empty() {
                set_mkl_blocked(format!(
                    "MKL redist incomplete (missing {}). Calling it would abort the process \
                     with no Axia error (Intel MKL FATAL ERROR / OpenMP). {inv}. \
                     Copy the full MKL bin redist next to axia.exe, or use --solver faer.",
                    missing.join(", ")
                ));
            } else {
                set_mkl_diag(format!(
                    "runtime {} → {} for wrapper; {inv}",
                    found.display(),
                    shim.display()
                ));
            }
            if found.file_name() != shim.file_name() {
                eprintln!(
                    "axia: MKL runtime {} (wrapper expects {expected}; using {})",
                    found.display(),
                    shim.display()
                );
            }
        }
        Err(e) => {
            set_mkl_blocked(format!(
                "found {} but could not create {expected}: {e}. {}",
                found.display(),
                mkl_dir_inventory(&dir)
            ));
        }
    }
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
    if csr.n == 0 {
        return Ok(Vec::new());
    }
    let n = csr.n as i32;
    let (a, ia, ja) = csr_upper_1based(csr)?;
    let try_type = |mtype: MatrixType| -> Result<Vec<f64>> {
        let inner = || -> Result<Vec<f64>> {
            let mut b = rhs.to_vec();
            let mut x = vec![0.0; csr.n];
            let mut ps = S::new().map_err(|e| FemError(e.to_string()))?;
            ps.set_matrix_type(mtype);
            ps.pardisoinit().map_err(|e| FemError(e.to_string()))?;
            // iparm[4]=0: do not read/write perm. Still allocate n slots —
            // pardiso-wrapper passes perm.as_mut_ptr() from an empty Vec.
            ps.set_perm(&vec![0i32; csr.n]);
            ps.set_iparm(4, 0);
            // iparm[26]=1: matrix checker (bad CSR → error code, not abort).
            ps.set_iparm(26, 1);
            ps.set_message_level(MessageLevel::Off);
            ps.set_phase(Phase::AnalysisNumFactSolveRefine);
            ps.pardiso(&a, &ia, &ja, &mut b, &mut x, n, 1)
                .map_err(|e| FemError(e.to_string()))?;
            Ok(x)
        };
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(inner)) {
            Ok(r) => r,
            Err(_) => err("PARDISO panicked (see axia: panic: … on stderr)"),
        }
    };
    match try_type(MatrixType::RealSymmetricPositiveDefinite) {
        Ok(x) => Ok(x),
        Err(_) => try_type(MatrixType::RealSymmetricIndefinite),
    }
}

fn probe_reexec_ok() -> bool {
    if std::env::var_os("AXIA_INTERNAL_MKL_PROBE").is_some() {
        return false;
    }
    if std::env::var_os("AXIA_SKIP_MKL_PROBE").is_some() {
        return false;
    }
    if cfg!(test) {
        return false;
    }
    let Ok(exe) = std::env::current_exe() else {
        return false;
    };
    let name = exe
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    name == "axia" || name.starts_with("axia-")
}

fn run_mkl_subprocess_probe() -> std::result::Result<(), String> {
    let exe = std::env::current_exe().map_err(|e| format!("current_exe: {e}"))?;
    let mut cmd = std::process::Command::new(&exe);
    cmd.arg("--internal-mkl-probe")
        .env("AXIA_INTERNAL_MKL_PROBE", "1")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        const CREATE_NO_WINDOW: u32 = 0x0800_0000;
        cmd.creation_flags(CREATE_NO_WINDOW);
    }
    let mut child = cmd
        .spawn()
        .map_err(|e| format!("could not spawn MKL self-test ({e})"))?;

    let mut stdout = child.stdout.take();
    let mut stderr = child.stderr.take();
    let t_out = std::thread::spawn(move || {
        let mut s = String::new();
        if let Some(ref mut r) = stdout {
            let _ = std::io::Read::read_to_string(r, &mut s);
        }
        s
    });
    let t_err = std::thread::spawn(move || {
        let mut s = String::new();
        if let Some(ref mut r) = stderr {
            let _ = std::io::Read::read_to_string(r, &mut s);
        }
        s
    });

    let start = std::time::Instant::now();
    let status = loop {
        match child.try_wait() {
            Ok(Some(st)) => break st,
            Ok(None) => {
                if start.elapsed() > std::time::Duration::from_secs(25) {
                    let _ = child.kill();
                    let _ = child.wait();
                    return Err(
                        "Intel MKL PARDISO self-test timed out (25s). \
                         Typical: missing CPU kernel DLL (mkl_avx2/mkl_def) or OpenMP deadlock."
                            .into(),
                    );
                }
                std::thread::sleep(std::time::Duration::from_millis(15));
            }
            Err(e) => return Err(format!("MKL self-test wait: {e}")),
        }
    };
    let out = t_out.join().unwrap_or_default();
    let err = t_err.join().unwrap_or_default();
    let log = format!("{out}{err}");
    let code = status.code();
    #[cfg(unix)]
    let signal = {
        use std::os::unix::process::ExitStatusExt;
        status.signal()
    };
    #[cfg(not(unix))]
    let signal: Option<i32> = None;

    if status.success() {
        return Ok(());
    }
    if code == Some(2) {
        let mut msg = if log.trim().is_empty() {
            "Intel MKL PARDISO is not available (self-test exit 2).".to_string()
        } else {
            log.trim().to_string()
        };
        if !msg.contains("axia:") {
            msg = format!("axia: {msg}");
        }
        return Err(msg);
    }
    let how = if let Some(s) = signal {
        format!("signal {s}")
    } else if let Some(c) = code {
        format!("exit {c}")
    } else {
        "aborted".to_string()
    };
    let mut msg = format!(
        "Intel MKL PARDISO crashed in a self-test ({how}) before the model was solved. \
         Typical causes: incomplete MKL redist (need mkl_core + mkl_avx2/mkl_def + threading DLL) \
         or an OpenMP conflict (libiomp5 vs libomp, e.g. Mecway). \
         Axia did not abort the job. Copy the full MKL bin folder next to axia.exe, \
         or use --solver faer."
    );
    if !log.trim().is_empty() {
        msg.push_str(" Output:\n");
        msg.push_str(log.trim());
    }
    Err(msg)
}

fn ensure_mkl_probe(complain: bool) -> bool {
    if !probe_reexec_ok() {
        return true;
    }
    static PROBE: std::sync::OnceLock<std::result::Result<(), String>> = std::sync::OnceLock::new();
    match PROBE.get_or_init(run_mkl_subprocess_probe) {
        Ok(()) => true,
        Err(msg) => {
            set_mkl_blocked(msg.clone());
            if complain {
                if msg.trim_start().starts_with("axia:") {
                    eprintln!("{msg}");
                } else {
                    eprintln!("axia: {msg}");
                }
            }
            false
        }
    }
}

fn try_mkl_inprocess(csr: &Csr, rhs: &[f64], complain: bool) -> Option<(Vec<f64>, String)> {
    #[cfg(target_arch = "x86_64")]
    {
        if pardiso_wrapper::MKLPardisoSolver::is_available() {
            announce("PARDISO (Intel MKL)");
            let _ = std::io::Write::flush(&mut std::io::stderr());
            match run_pardiso::<pardiso_wrapper::MKLPardisoSolver>(csr, rhs) {
                Ok(x) => return Some((x, "PARDISO (Intel MKL)".to_string())),
                Err(e) => {
                    if complain {
                        eprintln!("axia: PARDISO (Intel MKL) failed ({e})");
                    }
                    return None;
                }
            }
        } else if complain {
            use std::sync::atomic::{AtomicBool, Ordering};
            static HINT: AtomicBool = AtomicBool::new(false);
            if !HINT.swap(true, Ordering::Relaxed) {
                let extra = peek_mkl_diag().unwrap_or_default();
                let probe = find_mkl_runtime()
                    .or_else(|| {
                        exe_dir().and_then(|d| {
                            mkl_candidate_names()
                                .iter()
                                .map(|n| d.join(n))
                                .find(|p| p.is_file())
                        })
                    })
                    .map(|p| probe_load(&p))
                    .unwrap_or_default();
                eprintln!("axia: Intel MKL PARDISO did not load. {extra} {probe}");
            }
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = (csr, rhs);
        if complain {
            eprintln!("axia: Intel MKL PARDISO is only available on x86_64.");
        }
    }
    None
}

fn mkl_is_loaded() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        pardiso_wrapper::MKLPardisoSolver::is_available()
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

fn try_mkl(csr: &Csr, rhs: &[f64], complain: bool) -> Option<(Vec<f64>, String)> {
    prepare_mkl_env();
    if let Some(msg) = mkl_blocked_msg() {
        if complain {
            use std::sync::atomic::{AtomicBool, Ordering};
            static HINT: AtomicBool = AtomicBool::new(false);
            if !HINT.swap(true, Ordering::Relaxed) {
                eprintln!("axia: {msg}");
            }
        }
        return None;
    }
    if !mkl_is_loaded() {
        return try_mkl_inprocess(csr, rhs, complain);
    }
    if !ensure_mkl_probe(complain) {
        return None;
    }
    try_mkl_inprocess(csr, rhs, complain)
}

fn mkl_probe_csr() -> (Csr, Vec<f64>) {
    // 2×2 SPD diag(2,2) x = (2,2) → x = (1,1)
    let csr = Csr {
        n: 2,
        indptr: vec![0, 1, 2],
        indices: vec![0, 1],
        data: vec![2.0, 2.0],
    };
    (csr, vec![2.0, 2.0])
}

/// 0 = MKL solved a 2×2 system, non-zero = do not use MKL in the parent.
pub(crate) fn mkl_self_test() -> i32 {
    let (csr, b) = mkl_probe_csr();
    match try_mkl(&csr, &b, true) {
        Some((x, _)) if (x[0] - 1.0).abs() < 1e-6 && (x[1] - 1.0).abs() < 1e-6 => 0,
        Some((x, _)) => {
            eprintln!("axia: MKL probe produced unexpected x={x:?}");
            3
        }
        None => 2,
    }
}

fn try_panua(csr: &Csr, rhs: &[f64], complain: bool) -> Option<(Vec<f64>, String)> {
    if pardiso_wrapper::PanuaPardisoSolver::is_available() {
        let name = "PARDISO (Panua)";
        match run_pardiso::<pardiso_wrapper::PanuaPardisoSolver>(csr, rhs) {
            Ok(x) => return Some((x, name.to_string())),
            Err(e) => {
                if complain {
                    eprintln!("axia: {name} failed ({e})");
                }
            }
        }
    } else if complain {
        eprintln!("axia: Panua PARDISO is not on the loader path.");
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
        .map_err(|e| FemError(format!("Sparse-Matrixaufbau fehlgeschlagen ({e})")))
}

fn extract_col(x: &faer::Mat<f64>, n: usize) -> Vec<f64> {
    (0..n).map(|i| x[(i, 0)]).collect()
}

/// faer supernodal \(LL^\top\) (CHOLMOD-class, SPD). Pure Rust.
fn try_faer_llt(csr: &Csr, rhs: &[f64]) -> Result<Vec<f64>> {
    use faer::linalg::solvers::Solve;
    use faer::sparse::linalg::solvers::{Llt, SymbolicLlt};
    use faer::{Mat, Side};

    let mat = csr_to_faer(csr)?;
    let n = csr.n;
    let mut last = String::new();
    for side in [Side::Lower, Side::Upper] {
        match SymbolicLlt::<usize>::try_new(mat.symbolic(), side) {
            Ok(sym) => match Llt::try_new_with_symbolic(sym, mat.as_ref(), side) {
                Ok(llt) => {
                    let mut x = Mat::<f64>::from_fn(n, 1, |i, _| rhs[i]);
                    llt.solve_in_place(&mut x);
                    return Ok(extract_col(&x, n));
                }
                Err(e) => last = format!("{side:?} numeric: {e}"),
            },
            Err(e) => last = format!("{side:?} symbolic: {e}"),
        }
    }
    err(format!("faer LLT: {last}"))
}

/// faer supernodal \(LU\) for indefinite \(K\) (contact, some Newton steps).
fn try_faer_lu(csr: &Csr, rhs: &[f64]) -> Result<Vec<f64>> {
    use faer::linalg::solvers::Solve;
    use faer::sparse::linalg::solvers::{Lu, SymbolicLu};
    use faer::Mat;

    let mat = csr_to_faer(csr)?;
    let n = csr.n;
    let sym = SymbolicLu::<usize>::try_new(mat.symbolic())
        .map_err(|e| FemError(format!("faer LU (symbolic): {e}")))?;
    let lu = Lu::try_new_with_symbolic(sym, mat.as_ref())
        .map_err(|e| FemError(format!("faer LU: {e}")))?;
    let mut x = Mat::<f64>::from_fn(n, 1, |i, _| rhs[i]);
    lu.solve_in_place(&mut x);
    Ok(extract_col(&x, n))
}

fn try_faer(csr: &Csr, rhs: &[f64]) -> Result<(Vec<f64>, &'static str)> {
    match try_faer_llt(csr, rhs) {
        Ok(x) => Ok((x, "faer (supernodal LLT)")),
        Err(llt_e) => match try_faer_lu(csr, rhs) {
            Ok(x) => Ok((x, "faer (supernodal LU)")),
            Err(lu_e) => err(format!("{llt_e}; {lu_e}")),
        },
    }
}

fn pack(csr: &Csr, rhs: &[f64], x: Vec<f64>, name: &str) -> SparseSolve {
    let residual = residual(csr, &x, rhs);
    SparseSolve {
        x,
        name: name.to_string(),
        iters: 1,
        residual,
    }
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
    crate::backend::apply_env_solver();
    let want = crate::backend::sparse_backend();
    match want {
        SparseBackend::Cholesky | SparseBackend::Pcg => {
            return err("intern: dense/PCG laufen über linalg, nicht sparse_native");
        }
        SparseBackend::Mkl => {
            return match try_mkl(csr, rhs, true) {
                Some((x, name)) => {
                    announce(&name);
                    Ok(pack(csr, rhs, x, &name))
                }
                None => {
                    let extra = mkl_fail_detail();
                    if extra.is_empty() {
                        err("PARDISO (Intel MKL) nicht verfügbar oder Faktorisierung fehlgeschlagen.")
                    } else {
                        err(format!("PARDISO (Intel MKL) nicht verfügbar. {extra}"))
                    }
                }
            };
        }
        SparseBackend::Panua => {
            return match try_panua(csr, rhs, true) {
                Some((x, name)) => {
                    announce(&name);
                    Ok(pack(csr, rhs, x, &name))
                }
                None => err("PARDISO (Panua) nicht verfügbar oder Faktorisierung fehlgeschlagen."),
            };
        }
        SparseBackend::Pardiso => {
            if let Some((x, name)) = try_mkl(csr, rhs, true).or_else(|| try_panua(csr, rhs, true)) {
                announce(&name);
                return Ok(pack(csr, rhs, x, &name));
            }
            return err("PARDISO (MKL/Panua) nicht verfügbar.");
        }
        SparseBackend::Faer => {
            let (x, name) = try_faer(csr, rhs)?;
            announce(name);
            return Ok(pack(csr, rhs, x, name));
        }
        SparseBackend::Rivrs => {
            let name = "rivrs-sparse (LDLT)";
            announce(name);
            let x = try_rivrs(csr, rhs)?;
            return Ok(pack(csr, rhs, x, name));
        }
        SparseBackend::Auto => {}
    }

    if let Some((x, name)) = try_mkl(csr, rhs, true).or_else(|| try_panua(csr, rhs, false)) {
        announce(&name);
        return Ok(pack(csr, rhs, x, &name));
    }

    // SPD: faer supernodal LLT (PARDISO-class, no extra libs).
    // Indefinite (contact, some Newton): rivrs APTP LDLT, then faer LU.
    match try_faer_llt(csr, rhs) {
        Ok(x) => {
            announce("faer (supernodal LLT)");
            return Ok(pack(csr, rhs, x, "faer (supernodal LLT)"));
        }
        Err(_) => {}
    }

    match try_rivrs(csr, rhs) {
        Ok(x) => {
            announce("rivrs-sparse (LDLT)");
            return Ok(pack(csr, rhs, x, "rivrs-sparse (LDLT)"));
        }
        Err(e) => {
            eprintln!("axia: rivrs-sparse failed ({e}), trying faer LU");
        }
    }

    let x = try_faer_lu(csr, rhs)?;
    announce("faer (supernodal LU)");
    Ok(pack(csr, rhs, x, "faer (supernodal LU)"))
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
        assert!(!bin.join("libmkl_rt.dll").exists());
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn finds_dll_beside_exe_layout() {
        let tmp = std::env::temp_dir().join(format!(
            "axia-mkl-exe-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let dll = tmp.join("mkl_rt.dll");
        std::fs::write(&dll, b"fake-mkl").unwrap();
        let found = find_mkl_runtime_in(vec![tmp.clone()]).expect("dll next to exe");
        assert_eq!(found, dll);
        let shim = ensure_wrapper_named_library(&found).unwrap();
        assert_eq!(shim.file_name().unwrap(), wrapper_lib_name());
        assert_eq!(shim.parent().unwrap(), tmp.as_path());
        assert!(shim.is_file());
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn companion_name_detects_core_and_iomp() {
        assert!(is_mkl_companion_name("mkl_core.2.dll"));
        assert!(is_mkl_companion_name("mkl_intel_thread.2.dll"));
        assert!(is_mkl_companion_name("libiomp5md.dll"));
        assert!(is_mkl_companion_name("libmkl_rt.so.2"));
        assert!(is_mkl_companion_name("mkl_avx2.2.dll"));
        assert!(is_mkl_companion_name("mkl_def.2.dll"));
        assert!(!is_mkl_companion_name("axia.exe"));
        assert!(!is_mkl_companion_name("README.md"));
    }

    #[test]
    fn classifies_mkl_cpu_core_thread() {
        assert!(is_mkl_cpu_kernel_name("mkl_avx2.2.dll"));
        assert!(is_mkl_cpu_kernel_name("mkl_def.2.dll"));
        assert!(is_mkl_cpu_kernel_name("libmkl_avx512.so.2"));
        assert!(!is_mkl_cpu_kernel_name("mkl_core.2.dll"));
        assert!(is_mkl_core_name("mkl_core.2.dll"));
        assert!(is_mkl_core_name("libmkl_core.so.2"));
        assert!(is_mkl_thread_layer_name("mkl_sequential.2.dll"));
        assert!(is_mkl_thread_layer_name("mkl_intel_thread.2.dll"));
        assert!(is_iomp_name("libiomp5md.dll"));
        assert!(is_iomp_name("libiomp5.so"));
    }

    #[test]
    fn preflight_reports_missing_cpu_kernel() {
        let tmp = std::env::temp_dir().join(format!(
            "axia-mkl-preflight-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::fs::write(tmp.join("mkl_rt.dll"), b"fake").unwrap();
        let missing = missing_mkl_parts(&[tmp.clone()]);
        assert!(missing.iter().any(|s| s.contains("mkl_core")), "{missing:?}");
        assert!(missing.iter().any(|s| s.contains("avx2") || s.contains("def")), "{missing:?}");
        std::fs::write(tmp.join("mkl_core.2.dll"), b"fake").unwrap();
        std::fs::write(tmp.join("mkl_sequential.2.dll"), b"fake").unwrap();
        std::fs::write(tmp.join("mkl_avx2.2.dll"), b"fake").unwrap();
        let missing = missing_mkl_parts(&[tmp.clone()]);
        assert!(missing.is_empty(), "{missing:?}");
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn infer_mklroot_from_bin_dll() {
        let rt = std::path::PathBuf::from("/opt/intel/oneapi/mkl/latest/bin/mkl_rt.dll");
        let root = infer_mklroot(&rt).unwrap();
        assert_eq!(root, std::path::PathBuf::from("/opt/intel/oneapi/mkl/latest"));
        let rt = std::path::PathBuf::from("/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_rt.so.2");
        let root = infer_mklroot(&rt).unwrap();
        assert_eq!(root, std::path::PathBuf::from("/opt/intel/oneapi/mkl/latest"));
    }

    #[test]
    fn mkl_probe_matrix_is_spd_diag() {
        let (csr, b) = mkl_probe_csr();
        assert_eq!(csr.n, 2);
        assert_eq!(csr.indices, vec![0, 1]);
        assert_eq!(b, vec![2.0, 2.0]);
        let (a, ia, ja) = csr_upper_1based(&csr).unwrap();
        assert_eq!(ia, vec![1, 2, 3]);
        assert_eq!(ja, vec![1, 2]);
        assert_eq!(a, vec![2.0, 2.0]);
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
