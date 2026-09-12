//! Axia FEM command-line solver.
//!
//! CalculiX-style job runner: `axia job` reads `job.inp` and writes
//! `job.frd` + `job.dat`.

use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::process;

use axia_fem::{
    parse_model, parse_model_with_base, solve_native, solve_native_with_base, Model, SolveOutput,
};

const VERSION: &str = env!("CARGO_PKG_VERSION");

struct Args {
    job: Option<String>,
    input: Option<PathBuf>,
    output_stem: Option<PathBuf>,
    frd: Option<PathBuf>,
    dat: Option<PathBuf>,
    json: bool,
    check: bool,
    quiet: bool,
    verbose: bool,
    no_frd: bool,
    no_dat: bool,
    stdout_frd: bool,
}

fn print_help() {
    let exe = env::args().next().unwrap_or_else(|| "axia".into());
    println!(
        "\
Axia FEM {VERSION} — linear static solver (CalculiX INP / FRD / DAT)

Native sparse backends: PARDISO (MKL / Panua) if the library is on the
loader path, otherwise rivrs-sparse. The chosen solver is printed on stderr.

USAGE:
    {exe} [OPTIONS] <JOB>
    {exe} [OPTIONS] <JOB.inp>
    {exe} [OPTIONS] -              # read INP from stdin

ARGS:
    <JOB>                 Job name or path. Reads JOB.inp, writes JOB.frd and JOB.dat

OPTIONS:
    -i, --input <FILE>    Input INP (overrides JOB)
    -o, --output <STEM>   Output stem without extension (default: derived from input)
        --frd <FILE>      Write FRD to this path
        --dat <FILE>      Write DAT to this path
        --no-frd          Do not write a .frd file
        --no-dat          Do not write a .dat file
        --stdout          Write FRD to stdout (implies --quiet --no-dat)
        --json            Print solve statistics as JSON
        --check           Parse and report the mesh, do not solve
    -q, --quiet           Suppress the summary
    -v, --verbose         Extra diagnostics
    -h, --help            Show this help
    -V, --version         Show version

EXAMPLES:
    {exe} cantilever              # cantilever.inp → cantilever.frd / .dat
    {exe} model.inp
    {exe} -i deck.inp -o /tmp/run
    {exe} --check model.inp
    cat model.inp | {exe} - --json
"
    );
}

fn parse_args() -> Result<Args, String> {
    let mut args = Args {
        job: None,
        input: None,
        output_stem: None,
        frd: None,
        dat: None,
        json: false,
        check: false,
        quiet: false,
        verbose: false,
        no_frd: false,
        no_dat: false,
        stdout_frd: false,
    };
    let mut raw = env::args().skip(1);
    while let Some(a) = raw.next() {
        match a.as_str() {
            "-h" | "--help" => {
                print_help();
                process::exit(0);
            }
            "-V" | "--version" => {
                println!("axia {VERSION}");
                process::exit(0);
            }
            "-i" | "--input" => {
                args.input = Some(PathBuf::from(need(&mut raw, &a)?));
            }
            "-o" | "--output" => {
                args.output_stem = Some(PathBuf::from(need(&mut raw, &a)?));
            }
            "--frd" => args.frd = Some(PathBuf::from(need(&mut raw, &a)?)),
            "--dat" => args.dat = Some(PathBuf::from(need(&mut raw, &a)?)),
            "--json" => args.json = true,
            "--check" => args.check = true,
            "-q" | "--quiet" => args.quiet = true,
            "-v" | "--verbose" => args.verbose = true,
            "--no-frd" => args.no_frd = true,
            "--no-dat" => args.no_dat = true,
            "--stdout" => {
                args.stdout_frd = true;
                args.quiet = true;
                args.no_dat = true;
            }
            s if s.starts_with('-') => return Err(format!("unknown option '{s}'. Try --help.")),
            s => {
                if args.job.is_some() {
                    return Err("multiple job names given".into());
                }
                args.job = Some(s.to_string());
            }
        }
    }
    Ok(args)
}

fn need(it: &mut impl Iterator<Item = String>, flag: &str) -> Result<String, String> {
    it.next()
        .ok_or_else(|| format!("option {flag} requires a value"))
}

fn read_stdin() -> Result<String, String> {
    let mut s = String::new();
    io::stdin()
        .read_to_string(&mut s)
        .map_err(|e| format!("stdin: {e}"))?;
    if s.trim().is_empty() {
        return Err("stdin is empty".into());
    }
    Ok(s)
}

fn resolve_input(args: &Args) -> Result<(String, Option<PathBuf>), String> {
    if let Some(p) = &args.input {
        if p.as_os_str() == "-" {
            return read_stdin().map(|s| (s, None));
        }
        let text =
            fs::read_to_string(p).map_err(|e| format!("cannot read {}: {e}", p.display()))?;
        return Ok((text, Some(p.clone())));
    }
    match args.job.as_deref() {
        None => Err("missing job name. Try `axia --help`.".into()),
        Some("-") => read_stdin().map(|s| (s, None)),
        Some(job) => {
            let p = PathBuf::from(job);
            let path = if p.exists() {
                p
            } else {
                let with = PathBuf::from(format!("{job}.inp"));
                if with.exists() {
                    with
                } else {
                    return Err(format!("cannot find '{job}' or '{job}.inp'"));
                }
            };
            let text = fs::read_to_string(&path)
                .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
            Ok((text, Some(path)))
        }
    }
}

fn stem_from(args: &Args, input_path: Option<&Path>) -> PathBuf {
    if let Some(s) = &args.output_stem {
        return s.clone();
    }
    if let Some(p) = input_path {
        let mut s = p.to_path_buf();
        if s.extension().and_then(|e| e.to_str()).map(|e| e.eq_ignore_ascii_case("inp")) == Some(true)
        {
            s.set_extension("");
        }
        return s;
    }
    PathBuf::from("axia")
}

fn type_counts(model: &Model) -> BTreeMap<String, usize> {
    let mut m = BTreeMap::new();
    for e in &model.elements {
        *m.entry(e.kind.ccx_name().to_string()).or_insert(0) += 1;
    }
    m
}

fn format_types(model: &Model) -> String {
    let counts = type_counts(model);
    if counts.is_empty() {
        return "-".into();
    }
    counts
        .iter()
        .map(|(k, n)| format!("{n} {k}"))
        .collect::<Vec<_>>()
        .join(", ")
}

fn u_max(out: &SolveOutput) -> f64 {
    let mut umax = 0.0f64;
    for p in &out.u {
        umax = umax.max((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt());
    }
    umax
}

fn vm_range(out: &SolveOutput) -> (f64, f64) {
    let mut vmin = f64::MAX;
    let mut vmax = f64::MIN;
    for &m in &out.von_mises {
        vmin = vmin.min(m);
        vmax = vmax.max(m);
    }
    if vmin > vmax {
        (0.0, 0.0)
    } else {
        (vmin, vmax)
    }
}

fn print_model_summary(model: &Model, verbose: bool) {
    let heading = if model.heading.is_empty() {
        "-"
    } else {
        model.heading.as_str()
    };
    eprintln!("Axia FEM {VERSION}");
    eprintln!("  heading    {heading}");
    eprintln!("  nodes      {}", model.node_ids.len());
    eprintln!("  elements   {}  ({})", model.elements.len(), format_types(model));
    eprintln!("  dofs/node  {}", model.ndof_node());
    eprintln!("  dim        {}", model.dim);
    eprintln!("  bcs        {}", model.bcs.len());
    eprintln!("  cloads     {}", model.cloads.len());
    eprintln!("  dloads     {}", model.dloads.len());
    if verbose {
        for w in &model.warnings {
            eprintln!("  warning    {w}");
        }
    }
}

fn print_solve_summary(out: &SolveOutput) {
    let umax = u_max(out);
    let (vmin, vmax) = vm_range(out);
    eprintln!("  procedure  {}", out.procedure);
    eprintln!(
        "  dofs       {}  ({} free)",
        out.ndof, out.nfree
    );
    eprintln!("  solver     {}  ({} iter, r = {:.3e})", out.solver, out.iters, out.residual);
    eprintln!("  |u|_max    {umax:.6e}");
    eprintln!("  σ_vm       {vmin:.6e} … {vmax:.6e}");
    eprintln!("  time       {:.3} ms", out.time_ms);
}

fn stats_json(out: &SolveOutput) -> String {
    let umax = u_max(out);
    let (vmin, vmax) = vm_range(out);
    let types: serde_json::Map<String, serde_json::Value> = type_counts(&out.model)
        .into_iter()
        .map(|(k, n)| (k, serde_json::json!(n)))
        .collect();
    serde_json::json!({
        "ok": true,
        "version": VERSION,
        "heading": out.model.heading,
        "nnode": out.model.node_ids.len(),
        "nelem": out.model.elements.len(),
        "elements": types,
        "ndof": out.ndof,
        "nfree": out.nfree,
        "procedure": out.procedure,
        "ndofNode": out.model.ndof_node(),
        "solver": out.solver,
        "iterations": out.iters,
        "residual": out.residual,
        "timeMs": out.time_ms,
        "uMax": umax,
        "vmMin": vmin,
        "vmMax": vmax,
        "nbc": out.model.bcs.len(),
        "ncload": out.model.cloads.len(),
        "ndload": out.model.dloads.len(),
        "warnings": out.model.warnings,
    })
    .to_string()
}

fn write_file(path: &Path, contents: &str) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .map_err(|e| format!("cannot create {}: {e}", parent.display()))?;
        }
    }
    fs::write(path, contents).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    Ok(())
}

fn run() -> Result<(), String> {
    let args = parse_args()?;
    let (inp, input_path) = resolve_input(&args)?;
    let stem = stem_from(&args, input_path.as_deref());

    if args.check {
        let model = if let Some(p) = input_path.as_ref().and_then(|p| p.parent()) {
            parse_model_with_base(&inp, Some(p)).map_err(|e| e.to_string())?
        } else {
            parse_model(&inp).map_err(|e| e.to_string())?
        };
        if args.json {
            let types: serde_json::Map<String, serde_json::Value> = type_counts(&model)
                .into_iter()
                .map(|(k, n)| (k, serde_json::json!(n)))
                .collect();
            println!(
                "{}",
                serde_json::json!({
                    "ok": true,
                    "kind": "check",
                    "heading": model.heading,
                    "nnode": model.node_ids.len(),
                    "nelem": model.elements.len(),
                    "elements": types,
                    "ndofNode": model.ndof_node(),
                    "dim": model.dim,
                    "nbc": model.bcs.len(),
                    "ncload": model.cloads.len(),
                    "ndload": model.dloads.len(),
                    "warnings": model.warnings,
                })
            );
        } else if !args.quiet {
            print_model_summary(&model, args.verbose);
        }
        return Ok(());
    }

    let out = if let Some(p) = input_path.as_ref().and_then(|p| p.parent()) {
        solve_native_with_base(&inp, Some(p)).map_err(|e| e.to_string())?
    } else {
        solve_native(&inp).map_err(|e| e.to_string())?
    };

    if args.json {
        println!("{}", stats_json(&out));
    } else if !args.quiet {
        print_model_summary(&out.model, args.verbose);
        print_solve_summary(&out);
    }

    if args.stdout_frd {
        io::stdout()
            .write_all(out.frd.as_bytes())
            .map_err(|e| format!("stdout: {e}"))?;
        return Ok(());
    }

    if !args.no_frd {
        let path = args
            .frd
            .clone()
            .unwrap_or_else(|| stem.with_extension("frd"));
        write_file(&path, &out.frd)?;
        if !args.quiet && !args.json {
            eprintln!("  wrote      {}", path.display());
        }
    }
    if !args.no_dat {
        let path = args
            .dat
            .clone()
            .unwrap_or_else(|| stem.with_extension("dat"));
        write_file(&path, &out.dat)?;
        if !args.quiet && !args.json {
            eprintln!("  wrote      {}", path.display());
        }
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("axia: {e}");
        process::exit(1);
    }
}
