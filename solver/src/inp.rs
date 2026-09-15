use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::error::{err, Result};
use crate::model::{
    Amplitude, AnalysisStep, BeamSection, Boundary, Cflux, Cload, ContactPair, Coupling, Dflux,
    Dload, ElemKind, Element, Equation, Film, FluxKind, InitCond, InitKind, Material, Model,
    HyperKind, RigidBody, Surface, SurfaceInteraction, ThermalBc, Tie, Transform,
};

fn strip_comment(line: &str) -> &str {
    let t = line.trim();
    if t.starts_with("**") {
        return "";
    }
    line
}

fn parse_keyword(line: &str) -> (String, HashMap<String, String>) {
    let t = line.trim().trim_start_matches('\u{feff}');
    let upper = t.to_ascii_uppercase();
    let mut parts = upper.split(',');
    let kw = parts.next().unwrap_or("").trim().to_string();
    let mut params = HashMap::new();
    for p in parts {
        let p = p.trim();
        if p.is_empty() {
            continue;
        }
        if let Some((k, v)) = p.split_once('=') {
            params.insert(k.trim().to_string(), v.trim().to_string());
        } else {
            params.insert(p.to_string(), String::new());
        }
    }
    let kw = unglue_keyword(&kw);
    (kw, params)
}

/// CalculiX akzeptiert `*SOLIDSECTION` und `*SOLID SECTION`.
fn unglue_keyword(kw: &str) -> String {
    const TABLE: &[(&str, &str)] = &[
        ("*SOLIDSECTION", "*SOLID SECTION"),
        ("*SHELLSECTION", "*SHELL SECTION"),
        ("*MEMBRANESECTION", "*MEMBRANE SECTION"),
        ("*BEAMSECTION", "*BEAM SECTION"),
        ("*BEAMGENERALSECTION", "*BEAM GENERAL SECTION"),
        ("*RIGIDBODY", "*RIGID BODY"),
        ("*CONTACTPAIR", "*CONTACT PAIR"),
        ("*SURFACEINTERACTION", "*SURFACE INTERACTION"),
        ("*INITIALCONDITIONS", "*INITIAL CONDITIONS"),
        ("*ENDSTEP", "*END STEP"),
        ("*ELPRINT", "*EL PRINT"),
        ("*NODEPRINT", "*NODE PRINT"),
        ("*ELFILE", "*EL FILE"),
        ("*NODEFILE", "*NODE FILE"),
        ("*NODEOUTPUT", "*NODE OUTPUT"),
        ("*ELOUTPUT", "*EL OUTPUT"),
        ("*PHYSICALCONSTANTS", "*PHYSICAL CONSTANTS"),
        ("*TRANSFORM", "*TRANSFORM"),
    ];
    let u = kw.to_ascii_uppercase();
    for &(a, b) in TABLE {
        if u == a {
            return b.to_string();
        }
    }
    kw.to_string()
}

fn tokenize_data(line: &str) -> Vec<String> {
    line.split(|c: char| c == ',' || c.is_whitespace())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect()
}

fn parse_f64(s: &str) -> Result<f64> {
    parse_f64_inner(s).ok_or_else(|| crate::error::FemError(format!("Keine Zahl: {s}")))
}

/// CalculiX/Abaqus/Fortran floats: `1.d0`, `1.0D-3`, `1.`, `.5`, `1.0+3`.
fn parse_f64_inner(s: &str) -> Option<f64> {
    let raw = s.trim();
    if raw.is_empty() {
        return None;
    }
    let no_dot = raw.trim_end_matches('.');
    if let Ok(v) = no_dot.parse::<f64>() {
        return Some(v);
    }
    if let Ok(v) = raw.parse::<f64>() {
        return Some(v);
    }
    // Fortran D/d exponent → e
    let mut buf = String::with_capacity(raw.len());
    for c in raw.chars() {
        buf.push(if c == 'd' || c == 'D' { 'e' } else { c });
    }
    let b = buf.trim_end_matches('.');
    if let Ok(v) = b.parse::<f64>() {
        return Some(v);
    }
    // Fortran 1.0+3 / 1.0-3 (missing e) — only if a sign appears after a digit/dot.
    let bytes = b.as_bytes();
    let mut split = None;
    for i in 1..bytes.len() {
        if (bytes[i] == b'+' || bytes[i] == b'-') && (bytes[i - 1].is_ascii_digit() || bytes[i - 1] == b'.')
        {
            split = Some(i);
            break;
        }
    }
    if let Some(i) = split {
        let mut eform = String::with_capacity(b.len() + 1);
        eform.push_str(&b[..i]);
        eform.push('e');
        eform.push_str(&b[i..]);
        if let Ok(v) = eform.parse::<f64>() {
            return Some(v);
        }
    }
    None
}

fn parse_i32(s: &str) -> Result<i32> {
    let t = s.trim();
    if let Ok(v) = t.parse::<i32>() {
        return Ok(v);
    }
    if let Some(f) = parse_f64_inner(t) {
        return Ok(f.round() as i32);
    }
    err(format!("Keine ganze Zahl: {s}"))
}

fn looks_like_number(s: &str) -> bool {
    parse_f64_inner(s).is_some()
}

fn is_int_token(s: &str) -> bool {
    s.trim().parse::<i32>().is_ok()
}

/// `P`, `P2`, `P2NP`, `P3NU` → (face, nodal/nonuniform).
fn parse_dload_face(typ: &str) -> Option<(i32, bool)> {
    if typ == "P" {
        return Some((0, false));
    }
    if !typ.starts_with('P') {
        return None;
    }
    let rest = &typ[1..];
    let mut digits = String::new();
    for c in rest.chars() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else {
            break;
        }
    }
    let nodal = rest.contains("NP") || rest.contains("NU");
    if digits.is_empty() {
        return if nodal { Some((1, true)) } else { None };
    }
    let face: i32 = digits.parse().unwrap_or(0);
    Some((face, nodal))
}

fn ensure_node(model: &mut Model, id: i32) {
    if !model.node_ids.contains(&id) {
        model.node_ids.push(id);
        model.coords.push([0.0, 0.0, 0.0]);
        model.warn(format!("Knoten {id} angelegt (REF NODE ohne *NODE)."));
    }
}

fn is_keyword_line(line: &str) -> bool {
    let t = line.trim();
    t.starts_with('*') && !t.starts_with("**")
}

/// Collect data tokens from lines starting at `i` until the next keyword.
/// Returns (tokens, index of next keyword or end).
fn collect_tokens(lines: &[&str], mut i: usize) -> (Vec<String>, usize) {
    let mut tokens = Vec::new();
    while i < lines.len() {
        let raw = strip_comment(lines[i]);
        if raw.trim().is_empty() {
            i += 1;
            continue;
        }
        if is_keyword_line(raw) {
            break;
        }
        tokens.extend(tokenize_data(raw));
        i += 1;
    }
    (tokens, i)
}

fn next_nonempty(lines: &[&str], mut i: usize) -> usize {
    while i < lines.len() {
        let raw = strip_comment(lines[i]);
        if !raw.trim().is_empty() {
            return i;
        }
        i += 1;
    }
    i
}

pub fn parse(inp: &str) -> Result<Model> {
    parse_with_base(inp, None)
}

pub fn parse_with_base(inp: &str, base: Option<&Path>) -> Result<Model> {
    let expanded = expand_includes(inp, base, 0)?;
    parse_expanded(&expanded)
}

fn expand_includes(inp: &str, base: Option<&Path>, depth: usize) -> Result<String> {
    if depth > 16 {
        return err("*INCLUDE: Verschachtelung zu tief.");
    }
    let mut out = String::with_capacity(inp.len());
    for line in inp.lines() {
        let raw = strip_comment(line);
        if raw.trim().is_empty() {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if is_keyword_line(raw) {
            let (kw, params) = parse_keyword(raw);
            if kw == "*INCLUDE" {
                let file = include_filename(raw, line).ok_or_else(|| {
                    crate::error::FemError("*INCLUDE ohne INPUT=".into())
                })?;
                let path = resolve_include(base, &file);
                #[cfg(target_arch = "wasm32")]
                {
                    let _ = path;
                    out.push_str(&format!("** INCLUDE skipped in WASM: {file}\n"));
                    continue;
                }
                #[cfg(not(target_arch = "wasm32"))]
                {
                    let text = std::fs::read_to_string(&path).map_err(|e| {
                        crate::error::FemError(format!(
                            "*INCLUDE kann '{}' nicht lesen: {e}",
                            path.display()
                        ))
                    })?;
                    let nested_base = path.parent().map(|p| p.to_path_buf());
                    let nested = expand_includes(&text, nested_base.as_deref(), depth + 1)?;
                    out.push_str(&nested);
                    if !out.ends_with('\n') {
                        out.push('\n');
                    }
                    continue;
                }
            }
        }
        out.push_str(line);
        out.push('\n');
    }
    Ok(out)
}

fn resolve_include(base: Option<&Path>, file: &str) -> PathBuf {
    let p = PathBuf::from(file.trim().trim_matches('"').trim_matches('\''));
    if p.is_absolute() {
        p
    } else if let Some(b) = base {
        b.join(p)
    } else {
        p
    }
}

fn include_filename(raw_upper_line: &str, original: &str) -> Option<String> {
    let _ = raw_upper_line;
    let upper = original.to_ascii_uppercase();
    let key = if let Some(i) = upper.find("INPUT") {
        i
    } else {
        return None;
    };
    let rest = original[key + 5..].trim().trim_start_matches('=').trim();
    let token = rest.split([',', ' ', '\t']).find(|s| !s.is_empty())?;
    Some(token.trim_matches('"').trim_matches('\'').to_string())
}

fn parse_expanded(inp: &str) -> Result<Model> {
    let mut model = Model::new();
    let owned: Vec<String> = inp.lines().map(|l| l.to_string()).collect();
    let lines: Vec<&str> = owned.iter().map(|s| s.as_str()).collect();
    let n = lines.len();
    let mut i = 0;
    let mut current_material: Option<String> = None;
    let mut current_coupling: Option<usize> = None;
    let mut current_interaction: Option<String> = None;
    let mut saw_step = false;
    let mut step_open = false;
    let mut step_bc_from = 0usize;
    let mut step_cload_from = 0usize;
    let mut step_dload_from = 0usize;

    while i < n {
        let raw = strip_comment(lines[i]);
        if raw.trim().is_empty() {
            i += 1;
            continue;
        }
        if !is_keyword_line(raw) {
            i += 1;
            continue;
        }
        let (kw, params) = parse_keyword(raw);
        match kw.as_str() {
            "*HEADING" => {
                i += 1;
                i = next_nonempty(&lines, i);
                if i < n && !is_keyword_line(lines[i]) {
                    model.heading = lines[i].trim().to_string();
                    i += 1;
                }
            }
            "*NODE" => {
                let nset = params.get("NSET").cloned();
                i += 1;
                let mut added = Vec::new();
                // Per-line: ccx mixes `id`, `id,x,y` and `id,x,y,z` in one *NODE block.
                while i < n {
                    let raw = strip_comment(lines[i]);
                    if raw.trim().is_empty() {
                        i += 1;
                        continue;
                    }
                    if is_keyword_line(raw) {
                        break;
                    }
                    let toks = tokenize_data(raw);
                    i += 1;
                    if toks.is_empty() {
                        continue;
                    }
                    let id = parse_i32(&toks[0])?;
                    let x = if toks.len() > 1 {
                        parse_f64_inner(&toks[1]).unwrap_or(0.0)
                    } else {
                        0.0
                    };
                    let y = if toks.len() > 2 {
                        parse_f64_inner(&toks[2]).unwrap_or(0.0)
                    } else {
                        0.0
                    };
                    let z = if toks.len() > 3 {
                        parse_f64_inner(&toks[3]).unwrap_or(0.0)
                    } else {
                        0.0
                    };
                    model.node_ids.push(id);
                    model.coords.push([x, y, z]);
                    added.push(id);
                }
                if let Some(ns) = nset {
                    model.nsets.entry(ns).or_default().extend(added);
                }
            }
            "*ELEMENT" => {
                let typ = params.get("TYPE").cloned().ok_or_else(|| {
                    crate::error::FemError("*ELEMENT ohne TYPE=".into())
                })?;
                let kind = ElemKind::from_ccx(&typ).ok_or_else(|| {
                    crate::error::FemError(format!(
                        "Nicht unterstützter Elementtyp {typ}."
                    ))
                })?;
                let elset = params
                    .get("ELSET")
                    .cloned()
                    .unwrap_or_else(|| "EALL".into());
                let nn = kind.nnodes();
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut k = 0;
                let mut added = Vec::new();
                while k + nn < toks.len() {
                    let id = parse_i32(&toks[k])?;
                    k += 1;
                    let mut nodes = Vec::with_capacity(nn);
                    for _ in 0..nn {
                        nodes.push(parse_i32(&toks[k])?);
                        k += 1;
                    }
                    model.elements.push(Element {
                        id,
                        kind,
                        nodes,
                        elset: elset.clone(),
                    });
                    added.push(id);
                }
                if k != toks.len() {
                    model.warn(format!(
                        "ELEMENT {typ}: {} Zahlen übrig, ignoriert.",
                        toks.len() - k
                    ));
                }
                model.elsets.entry(elset).or_default().extend(added);
            }
            "*NSET" => {
                let name = params.get("NSET").cloned().ok_or_else(|| {
                    crate::error::FemError("*NSET ohne NSET=".into())
                })?;
                let generate = params.contains_key("GENERATE");
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut ids = Vec::new();
                if generate {
                    let mut k = 0;
                    while k + 1 < toks.len() {
                        let a = parse_i32(&toks[k])?;
                        let b = parse_i32(&toks[k + 1])?;
                        let step = if k + 2 < toks.len() {
                            let s = parse_i32(&toks[k + 2]).unwrap_or(1);
                            k += 3;
                            if s == 0 { 1 } else { s }
                        } else {
                            k += 2;
                            1
                        };
                        let mut v = a;
                        if step > 0 {
                            while v <= b {
                                ids.push(v);
                                v += step;
                            }
                        } else {
                            while v >= b {
                                ids.push(v);
                                v += step;
                            }
                        }
                    }
                } else {
                    for t in &toks {
                        match parse_i32(t) {
                            Ok(id) => ids.push(id),
                            Err(_) => model.warnings.push(format!(
                                "__NS__|{name}|{}",
                                t.to_ascii_uppercase()
                            )),
                        }
                    }
                }
                model.nsets.entry(name).or_default().extend(ids);
            }
            "*ELSET" => {
                let name = params.get("ELSET").cloned().ok_or_else(|| {
                    crate::error::FemError("*ELSET ohne ELSET=".into())
                })?;
                let generate = params.contains_key("GENERATE");
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut ids = Vec::new();
                if generate {
                    let mut k = 0;
                    while k + 1 < toks.len() {
                        let a = parse_i32(&toks[k])?;
                        let b = parse_i32(&toks[k + 1])?;
                        let step = if k + 2 < toks.len() {
                            let s = parse_i32(&toks[k + 2]).unwrap_or(1);
                            k += 3;
                            if s == 0 { 1 } else { s }
                        } else {
                            k += 2;
                            1
                        };
                        let mut v = a;
                        if step > 0 {
                            while v <= b {
                                ids.push(v);
                                v += step;
                            }
                        } else {
                            while v >= b {
                                ids.push(v);
                                v += step;
                            }
                        }
                    }
                } else {
                    for t in &toks {
                        match parse_i32(t) {
                            Ok(id) => ids.push(id),
                            Err(_) => model.warnings.push(format!(
                                "__ES__|{name}|{}",
                                t.to_ascii_uppercase()
                            )),
                        }
                    }
                }
                model.elsets.entry(name).or_default().extend(ids);
            }
            "*MATERIAL" => {
                let name = params
                    .get("NAME")
                    .cloned()
                    .unwrap_or_else(|| "MATERIAL-1".into());
                model.materials.entry(name.clone()).or_default();
                current_material = Some(name);
                i += 1;
            }
            "*ELASTIC" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if toks.len() < 2 {
                    return err("*ELASTIC benötigt E und nu");
                }
                let typ = params
                    .get("TYPE")
                    .map(|s| s.to_ascii_uppercase())
                    .unwrap_or_default();
                let (e, nu) = if typ.contains("ORTHO") || typ.contains("ANISO") {
                    // ccx TYPE=ORTHO/ANISO: D11,D12,... possibly temperature tables
                    // starting with a zero row. Skip leading zeros; isotropic approx.
                    let vals: Vec<f64> = toks.iter().filter_map(|t| parse_f64_inner(t)).collect();
                    let start = vals.iter().position(|v| v.abs() > 1e-18).unwrap_or(0);
                    let d11 = vals.get(start).copied().unwrap_or(0.0);
                    let d12 = vals.get(start + 1).copied().unwrap_or(0.0);
                    let e = if d11.abs() > 1e-18 { d11 } else { 1.0 };
                    let nu = if e.abs() > 1e-18 {
                        (d12 / e).clamp(-0.49, 0.49)
                    } else {
                        0.3
                    };
                    (e, nu)
                } else if typ.contains("ENGINEERING") && toks.len() >= 4 {
                    let e = parse_f64(&toks[0])?;
                    let nu = parse_f64(&toks[3])?;
                    let e = if e.abs() > 1e-18 {
                        e
                    } else {
                        toks.iter()
                            .filter_map(|t| parse_f64_inner(t))
                            .find(|v| v.abs() > 1e-18)
                            .unwrap_or(1.0)
                    };
                    (e, nu.clamp(-0.49, 0.49))
                } else {
                    let mut e = parse_f64(&toks[0])?;
                    let mut nu = parse_f64(&toks[1])?;
                    if e.abs() <= 1e-18 {
                        if let Some(v) = toks
                            .iter()
                            .filter_map(|t| parse_f64_inner(t))
                            .find(|v| v.abs() > 1e-18)
                        {
                            e = v;
                            model.warn("*ELASTIC: erste Zeile E=0, temperaturabhängige Zeile verwendet.");
                        } else {
                            e = 1.0;
                            model.warn("*ELASTIC: E=0 — Dummy E=1 verwendet.");
                        }
                    }
                    if !(-0.49..=0.49).contains(&nu) {
                        nu = nu.clamp(-0.49, 0.49);
                    }
                    (e, nu)
                };
                let name = current_material
                    .clone()
                    .unwrap_or_else(|| "MATERIAL-1".into());
                let m = model.materials.entry(name).or_default();
                m.e = e;
                m.nu = nu;
                if typ.contains("ORTHO") || typ.contains("ANISO") || typ.contains("ENGINEERING") {
                    model.warn(format!(
                        "*ELASTIC, TYPE={typ} als isotrop angenähert (E={e}, nu={nu})."
                    ));
                }
            }
            "*DENSITY" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if toks.is_empty() {
                    return err("*DENSITY ohne Wert");
                }
                let rho = parse_f64(&toks[0])?;
                let name = current_material
                    .clone()
                    .unwrap_or_else(|| "MATERIAL-1".into());
                model.materials.entry(name).or_default().density = rho;
            }
            "*SOLID SECTION" | "*SHELL SECTION" | "*MEMBRANE SECTION" => {
                let elset = params.get("ELSET").cloned().unwrap_or_else(|| "EALL".into());
                if let Some(mat) = params.get("MATERIAL") {
                    model.elset_material.insert(elset.clone(), mat.clone());
                }
                let composite = params.contains_key("COMPOSITE");
                if composite {
                    i += 1;
                    let mut tsum = 0.0;
                    let mut mat_name: Option<String> = None;
                    while i < n {
                        let raw = strip_comment(lines[i]);
                        if raw.trim().is_empty() {
                            i += 1;
                            continue;
                        }
                        if is_keyword_line(raw) {
                            break;
                        }
                        let parts: Vec<&str> = raw.split(',').map(|s| s.trim()).collect();
                        if let Some(th) = parts.first().and_then(|s| parse_f64_inner(s)) {
                            tsum += th.abs();
                        }
                        if let Some(m) = parts.get(2).filter(|s| !s.is_empty()) {
                            if mat_name.is_none() {
                                mat_name = Some(m.to_ascii_uppercase());
                            }
                        }
                        i += 1;
                    }
                    if let Some(m) = mat_name {
                        model.elset_material.entry(elset.clone()).or_insert(m);
                    }
                    model.elset_thickness.insert(elset, tsum.max(1e-12));
                } else {
                    let (toks, ni) = collect_tokens(&lines, i + 1);
                    i = ni;
                    let t = if !toks.is_empty() {
                        parse_f64(&toks[0]).unwrap_or(1.0)
                    } else {
                        1.0
                    };
                    model.elset_thickness.insert(elset, t);
                }
            }
            "*BEAM SECTION" | "*BEAM GENERAL SECTION" => {
                let elset = params.get("ELSET").cloned().unwrap_or_else(|| "EALL".into());
                if let Some(mat) = params.get("MATERIAL") {
                    model.elset_material.insert(elset.clone(), mat.clone());
                }
                let sectyp = params
                    .get("SECTION")
                    .cloned()
                    .unwrap_or_else(|| "RECT".into());
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if toks.is_empty() {
                    return err("*BEAM SECTION ohne Querschnittswerte");
                }
                let sec = parse_beam_section(&sectyp, &toks)?;
                model.elset_beam.insert(elset, sec);
            }
            "*BOUNDARY" => {
                if op_is_new(&params) {
                    step_bc_from = model.bcs.len();
                }
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut k = 0;
                while k < toks.len() {
                    let name = toks[k].to_ascii_uppercase();
                    k += 1;
                    if k >= toks.len() {
                        break;
                    }
                    let d1 = parse_i32(&toks[k])? as usize;
                    k += 1;
                    let d2 = if k < toks.len() && parse_i32(&toks[k]).is_ok() {
                        let v = parse_i32(&toks[k])? as usize;
                        k += 1;
                        v
                    } else {
                        d1
                    };
                    let first = d1.max(1);
                    let last = d2.max(first);
                    let thermal = first == 11 || last == 11;
                    let val = if k < toks.len() {
                        let peek = &toks[k];
                        if thermal && looks_like_number(peek) {
                            let v = parse_f64_inner(peek).unwrap_or(0.0);
                            k += 1;
                            v
                        } else if looks_like_number(peek) {
                            let is_pure_int = is_int_token(peek);
                            if is_pure_int && k + 1 < toks.len() {
                                let maybe_dof = parse_i32(&toks[k + 1]).ok();
                                if maybe_dof == Some(1)
                                    || maybe_dof == Some(2)
                                    || maybe_dof == Some(3)
                                    || maybe_dof == Some(4)
                                    || maybe_dof == Some(5)
                                    || maybe_dof == Some(6)
                                    || maybe_dof == Some(11)
                                {
                                    0.0
                                } else {
                                    let v = parse_f64_inner(peek).unwrap_or(0.0);
                                    k += 1;
                                    v
                                }
                            } else {
                                let v = parse_f64_inner(peek).unwrap_or(0.0);
                                k += 1;
                                v
                            }
                        } else {
                            0.0
                        }
                    } else {
                        0.0
                    };
                    let nodes = if let Ok(id) = parse_i32(&name) {
                        vec![id]
                    } else {
                        // nsets not yet complete — store as synthetic later
                        model
                            .nsets
                            .get(&name)
                            .cloned()
                            .unwrap_or_else(|| vec![])
                            .into_iter()
                            .chain(if name.chars().all(|c| c.is_ascii_digit() || c == '-') {
                                vec![]
                            } else {
                                Vec::new()
                            })
                            .collect()
                    };
                    if nodes.is_empty() && parse_i32(&name).is_err() {
                    }
                    if parse_i32(&name).is_ok() {
                        let id = parse_i32(&name)?;
                        for d in first..=last {
                            if d >= 1 && d <= 6 {
                                model.bcs.push(Boundary {
                                    node: id,
                                    dof: d - 1,
                                    value: val,
                                });
                            } else if d == 11 {
                                model.thermal_bcs.push(ThermalBc { node: id, value: val });
                            }
                        }
                    } else if let Some(nodes) = model.nsets.get(&name).cloned().filter(|v| !v.is_empty()) {
                        for id in nodes {
                            for d in first..=last {
                                if d >= 1 && d <= 6 {
                                    model.bcs.push(Boundary {
                                        node: id,
                                        dof: d - 1,
                                        value: val,
                                    });
                                } else if d == 11 {
                                    model.thermal_bcs.push(ThermalBc { node: id, value: val });
                                }
                            }
                        }
                    } else {
                        for d in first..=last {
                            if d >= 1 && d <= 6 {
                                model.warnings.push(format!(
                                    "__BC__|{name}|{}|{}|{val}",
                                    d - 1,
                                    d - 1
                                ));
                            } else if d == 11 {
                                model.warnings.push(format!("__TB__|{name}|{val}"));
                            }
                        }
                    }
                }
            }
            "*CLOAD" => {
                if op_is_new(&params) {
                    step_cload_from = model.cloads.len();
                }
                let amp_name = params
                    .get("AMPLITUDE")
                    .cloned()
                    .unwrap_or_default()
                    .to_ascii_uppercase();
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut k = 0;
                while k + 2 < toks.len() {
                    let name = toks[k].to_ascii_uppercase();
                    let dof = parse_i32(&toks[k + 1])? as usize;
                    let mag = parse_f64(&toks[k + 2])?;
                    k += 3;
                    if dof < 1 || dof > 6 {
                        model.warn(format!("CLOAD Freiheitsgrad {dof} ignoriert"));
                        continue;
                    }
                    if let Ok(id) = parse_i32(&name) {
                        model.cloads.push(Cload {
                            node: id,
                            dof: dof - 1,
                            mag,
                            amplitude: amp_name.clone(),
                        });
                    } else if let Some(nodes) = model.nsets.get(&name).cloned().filter(|v| !v.is_empty()) {
                        for id in nodes {
                            model.cloads.push(Cload {
                                node: id,
                                dof: dof - 1,
                                mag,
                                amplitude: amp_name.clone(),
                            });
                        }
                    } else {
                        model.warnings.push(format!(
                            "__CL__|{name}|{}|{mag}|{amp_name}",
                            dof - 1
                        ));
                    }
                }
            }
            "*DLOAD" => {
                if op_is_new(&params) {
                    step_dload_from = model.dloads.len();
                }
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut k = 0;
                while k + 1 < toks.len() {
                    let name = toks[k].to_ascii_uppercase();
                    let typ = toks[k + 1].to_ascii_uppercase();
                    k += 2;
                    if typ == "GRAV" {
                        if k + 3 >= toks.len() {
                            return err("*DLOAD GRAV benötigt Betrag und Richtung");
                        }
                        let mag = parse_f64(&toks[k])?;
                        let dx = parse_f64(&toks[k + 1])?;
                        let dy = parse_f64(&toks[k + 2])?;
                        let dz = parse_f64(&toks[k + 3])?;
                        k += 4;
                        model.dloads.push(Dload::Grav {
                            mag,
                            dir: [dx, dy, dz],
                        });
                    } else if typ == "CENTRIF" {
                        if k + 6 >= toks.len() {
                            return err(
                                "*DLOAD CENTRIF braucht ω² und zwei Achspunkte (7 Zahlen).",
                            );
                        }
                        let omega2 = parse_f64(&toks[k])?;
                        let p1 = [
                            parse_f64(&toks[k + 1])?,
                            parse_f64(&toks[k + 2])?,
                            parse_f64(&toks[k + 3])?,
                        ];
                        let p2 = [
                            parse_f64(&toks[k + 4])?,
                            parse_f64(&toks[k + 5])?,
                            parse_f64(&toks[k + 6])?,
                        ];
                        k += 7;
                        if parse_i32(&name).is_ok() {
                            let id = parse_i32(&name)?;
                            model.dloads.push(Dload::Centrif {
                                omega2,
                                p1,
                                p2,
                                elems: vec![id],
                            });
                        } else if let Some(elems) =
                            model.elsets.get(&name).cloned().filter(|v| !v.is_empty())
                        {
                            model.dloads.push(Dload::Centrif {
                                omega2,
                                p1,
                                p2,
                                elems,
                            });
                        } else {
                            model.warnings.push(format!(
                                "__CFUG__|{name}|{omega2}|{}|{}|{}|{}|{}|{}",
                                p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]
                            ));
                        }
                    } else if typ == "PX" || typ == "PY" || typ == "PZ" {
                        if k >= toks.len() {
                            return err("*DLOAD PX/PY/PZ ohne Betrag");
                        }
                        let mag = parse_f64(&toks[k])?;
                        k += 1;
                        let dir = match typ.as_str() {
                            "PX" => [1.0, 0.0, 0.0],
                            "PY" => [0.0, 1.0, 0.0],
                            _ => [0.0, 0.0, 1.0],
                        };
                        if let Ok(id) = parse_i32(&name) {
                            model.dloads.push(Dload::BeamGlobal { elem: id, dir, mag });
                        } else {
                            model.warnings.push(format!(
                                "__BG__|{name}|{}|{}|{}|{mag}",
                                dir[0], dir[1], dir[2]
                            ));
                        }
                    } else if let Some((face, nodal)) = parse_dload_face(&typ) {
                        if k >= toks.len() {
                            return err("*DLOAD P ohne Betrag");
                        }
                        let mut mag;
                        let mut ref_node: Option<i32> = None;
                        if nodal && is_int_token(&toks[k]) {
                            ref_node = parse_i32(&toks[k]).ok();
                            mag = 1.0;
                            k += 1;
                            let mut vals = Vec::new();
                            while k < toks.len() && looks_like_number(&toks[k]) && !is_int_token(&toks[k])
                            {
                                vals.push(parse_f64_inner(&toks[k]).unwrap_or(0.0));
                                k += 1;
                            }
                            if !vals.is_empty() {
                                mag = vals.iter().sum::<f64>() / vals.len() as f64;
                                ref_node = None;
                            }
                        } else {
                            mag = parse_f64(&toks[k])?;
                            k += 1;
                        }
                        if let Some(nid) = ref_node {
                            if parse_i32(&name).is_ok() {
                                model.warnings.push(format!("__DLN__|{name}|{face}|{nid}"));
                            } else {
                                model.warnings.push(format!("__DLN__|{name}|{face}|{nid}"));
                            }
                        } else if let Ok(id) = parse_i32(&name) {
                            model.dloads.push(Dload::Pressure {
                                elem: id,
                                face,
                                mag,
                            });
                        } else if let Some(elems) =
                            model.elsets.get(&name).cloned().filter(|v| !v.is_empty())
                        {
                            for elem in elems {
                                model.dloads.push(Dload::Pressure { elem, face, mag });
                            }
                        } else {
                            model.warnings.push(format!("__DL__|{name}|{face}|{mag}"));
                        }
                    } else {
                        model.warn(format!("DLOAD-Typ {typ} nicht unterstützt"));
                        // skip one magnitude if present
                        if k < toks.len() && parse_f64(&toks[k]).is_ok() {
                            k += 1;
                        }
                    }
                }
            }
            "*STEP" => {
                if step_open {
                    push_step(
                        &mut model,
                        step_bc_from,
                        step_cload_from,
                        step_dload_from,
                    );
                }
                saw_step = true;
                step_open = true;
                if params.contains_key("NLGEOM") {
                    model.procedure = crate::model::Procedure::Static {
                        nlgeom: true,
                        increments: 1,
                        riks: false,
                    };
                }
                if let Some(inc) = params.get("INC") {
                    if let Ok(v) = parse_i32(inc) {
                        if v > 0 {
                            model.max_inc = v as usize;
                        }
                    }
                }
                i += 1;
            }
            "*STATIC" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                // CalculiX: initial Δt, period, min Δt, max Δt [, max increments for RIKS]
                let dt0 = toks
                    .get(0)
                    .and_then(|t| parse_f64(t).ok())
                    .filter(|v| *v > 0.0)
                    .unwrap_or(1.0);
                let period = toks
                    .get(1)
                    .and_then(|t| parse_f64(t).ok())
                    .filter(|v| *v > 0.0)
                    .unwrap_or(1.0);
                model.static_dt = dt0;
                model.static_period = period;
                let n_from_time = ((period / dt0).round() as usize).max(1);
                let inc = n_from_time.min(model.max_inc.max(1));
                let nlgeom = matches!(
                    model.procedure,
                    crate::model::Procedure::Static { nlgeom: true, .. }
                );
                let riks = params.contains_key("RIKS");
                if riks {
                    let mut ctrl = crate::model::RiksCtrl::default();
                    if let Some(v) = toks.get(0).and_then(|t| parse_f64(t).ok()) {
                        if v > 0.0 {
                            ctrl.dlam = v;
                        }
                    }
                    if let Some(v) = toks.get(1).and_then(|t| parse_f64(t).ok()) {
                        if v > 0.0 {
                            ctrl.period = v;
                        }
                    }
                    if let Some(v) = toks.get(2).and_then(|t| parse_f64(t).ok()) {
                        if v > 0.0 {
                            ctrl.dlam_min = v;
                        }
                    }
                    if let Some(v) = toks.get(3).and_then(|t| parse_f64(t).ok()) {
                        if v > 0.0 {
                            ctrl.dlam_max = v;
                        }
                    }
                    if let Some(v) = toks.get(4).and_then(|t| parse_i32(t).ok()) {
                        if v > 0 {
                            ctrl.max_inc = v as usize;
                        }
                    }
                    if ctrl.dlam_max < ctrl.dlam_min {
                        ctrl.dlam_max = ctrl.dlam_min;
                    }
                    if ctrl.dlam > ctrl.dlam_max {
                        ctrl.dlam_max = ctrl.dlam;
                    }
                    model.riks = Some(ctrl);
                    if !nlgeom {
                        model.warn("*STATIC, RIKS: NLGEOM wird eingeschaltet.");
                    }
                }
                model.procedure = crate::model::Procedure::Static {
                    nlgeom: nlgeom || riks,
                    increments: if riks {
                        model
                            .riks
                            .map(|c| c.max_inc)
                            .unwrap_or(inc)
                            .max(1)
                    } else {
                        inc.max(1)
                    },
                    riks,
                };
            }
            "*FREQUENCY" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let n = if !toks.is_empty() {
                    parse_i32(&toks[0]).unwrap_or(10).max(1) as usize
                } else {
                    10
                };
                model.procedure = crate::model::Procedure::Frequency { nmodes: n };
            }
            "*BUCKLE" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let n = if !toks.is_empty() {
                    parse_i32(&toks[0]).unwrap_or(1).max(1) as usize
                } else {
                    1
                };
                model.procedure = crate::model::Procedure::Buckle { nmodes: n };
            }
            "*HEAT TRANSFER" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let steady = params.contains_key("STEADY STATE")
                    || params.contains_key("STEADYSTATE");
                let dt = if !toks.is_empty() {
                    parse_f64(&toks[0]).unwrap_or(0.0)
                } else {
                    0.0
                };
                let period = if toks.len() >= 2 {
                    parse_f64(&toks[1]).unwrap_or(0.0)
                } else {
                    0.0
                };
                model.procedure = crate::model::Procedure::HeatTransfer {
                    steady: steady || dt <= 0.0 || period <= 0.0,
                    dt: dt.max(0.0),
                    period: period.max(0.0),
                };
                model.output_nt = true;
            }
            "*DYNAMIC" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let (dt, period) = match toks.len() {
                    0 => (0.01, 1.0),
                    1 => {
                        let p = parse_f64(&toks[0])?.abs().max(1e-16);
                        ((p / 20.0).max(1e-16), p)
                    }
                    _ => {
                        let mut dt = parse_f64(&toks[0]).unwrap_or(0.0);
                        let mut period = parse_f64(&toks[1]).unwrap_or(0.0);
                        if period <= 0.0 && dt > 0.0 {
                            period = dt;
                            dt = (period / 20.0).max(1e-16);
                        }
                        if dt <= 0.0 && period > 0.0 {
                            dt = (period / 20.0).max(1e-16);
                        }
                        if dt <= 0.0 || period <= 0.0 {
                            (0.01, 1.0)
                        } else {
                            (dt, period)
                        }
                    }
                };
                model.procedure = crate::model::Procedure::Dynamic { dt, period };
            }
            "*DAMPING" => {
                if let Some(a) = params.get("ALPHA").or_else(|| params.get("ALPHA")) {
                    if let Ok(v) = parse_f64(a) {
                        model.damp_alpha = v;
                    }
                }
                if let Some(b) = params.get("BETA") {
                    if let Ok(v) = parse_f64(b) {
                        model.damp_beta = v;
                    }
                }
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if model.damp_alpha == 0.0 && !toks.is_empty() {
                    model.damp_alpha = parse_f64(&toks[0]).unwrap_or(0.0);
                }
                if model.damp_beta == 0.0 && toks.len() >= 2 {
                    model.damp_beta = parse_f64(&toks[1]).unwrap_or(0.0);
                }
            }
            "*CONDUCTIVITY" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if toks.is_empty() {
                    return err("*CONDUCTIVITY ohne Wert");
                }
                let k = parse_f64(&toks[0])?;
                let name = current_material
                    .clone()
                    .unwrap_or_else(|| "MATERIAL-1".into());
                model.materials.entry(name).or_default().conductivity = k;
            }
            "*SPECIFIC HEAT" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if toks.is_empty() {
                    return err("*SPECIFIC HEAT ohne Wert");
                }
                let cp = parse_f64(&toks[0])?;
                let name = current_material
                    .clone()
                    .unwrap_or_else(|| "MATERIAL-1".into());
                model.materials.entry(name).or_default().specific_heat = cp;
            }
            "*CFLUX" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut k = 0;
                while k + 1 < toks.len() {
                    let name = toks[k].to_ascii_uppercase();
                    k += 1;
                    // optional dof 11
                    if k < toks.len() {
                        if let Ok(d) = parse_i32(&toks[k]) {
                            if d == 11 || d == 0 {
                                k += 1;
                            }
                        }
                    }
                    if k >= toks.len() {
                        break;
                    }
                    let mag = parse_f64(&toks[k])?;
                    k += 1;
                    if let Ok(id) = parse_i32(&name) {
                        model.cfluxes.push(Cflux { node: id, mag });
                    } else {
                        model.warnings.push(format!("__CF__|{name}|{mag}"));
                    }
                }
            }
            "*DFLUX" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut k = 0;
                while k + 1 < toks.len() {
                    let name = toks[k].to_ascii_uppercase();
                    let typ = toks[k + 1].to_ascii_uppercase();
                    k += 2;
                    if k >= toks.len() {
                        return err("*DFLUX ohne Betrag");
                    }
                    let mag = parse_f64(&toks[k])?;
                    k += 1;
                    let kind = if typ == "BF" || typ == "BFNU" {
                        FluxKind::Body
                    } else if typ == "S" {
                        FluxKind::Face(0)
                    } else if typ.starts_with('S') {
                        let f = typ.trim_start_matches('S').parse::<i32>().unwrap_or(1);
                        FluxKind::Face(f)
                    } else {
                        model.warn(format!("DFLUX-Typ {typ} nicht unterstützt"));
                        continue;
                    };
                    if let Ok(id) = parse_i32(&name) {
                        model.dfluxes.push(Dflux {
                            elem: id,
                            kind,
                            mag,
                        });
                    } else {
                        let tag = match kind {
                            FluxKind::Body => 0,
                            FluxKind::Face(f) => f,
                        };
                        model
                            .warnings
                            .push(format!("__DF__|{name}|{tag}|{mag}"));
                    }
                }
            }
            "*FILM" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut k = 0;
                while k + 2 < toks.len() {
                    let name = toks[k].to_ascii_uppercase();
                    let typ = toks[k + 1].to_ascii_uppercase();
                    k += 2;
                    if k + 1 >= toks.len() {
                        return err("*FILM erwartet T∞ und h");
                    }
                    let tinf = parse_f64(&toks[k])?;
                    let h = parse_f64(&toks[k + 1])?;
                    k += 2;
                    let face = if typ == "F" {
                        0
                    } else if typ.starts_with('F') {
                        typ.trim_start_matches('F').parse::<i32>().unwrap_or(1)
                    } else {
                        1
                    };
                    if let Ok(id) = parse_i32(&name) {
                        model.films.push(Film {
                            elem: id,
                            face,
                            t_inf: tinf,
                            h,
                        });
                    } else {
                        model
                            .warnings
                            .push(format!("__FM__|{name}|{face}|{tinf}|{h}"));
                    }
                }
            }
            "*AMPLITUDE" => {
                let name = params
                    .get("NAME")
                    .cloned()
                    .unwrap_or_else(|| format!("A{}", model.amplitudes.len() + 1))
                    .to_ascii_uppercase();
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut pts = Vec::new();
                let mut k = 0;
                while k + 1 < toks.len() {
                    let t = parse_f64(&toks[k])?;
                    let v = parse_f64(&toks[k + 1])?;
                    k += 2;
                    pts.push((t, v));
                }
                model.amplitudes.push(Amplitude { name, points: pts });
            }
            "*INITIAL CONDITIONS" => {
                let ty = params
                    .get("TYPE")
                    .cloned()
                    .unwrap_or_else(|| "TEMPERATURE".into());
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let kind = match ty.as_str() {
                    "VELOCITY" => InitKind::Velocity,
                    "DISPLACEMENT" => InitKind::Displacement,
                    _ => InitKind::Temperature,
                };
                let mut k = 0;
                while k < toks.len() {
                    let name = toks[k].to_ascii_uppercase();
                    k += 1;
                    if kind == InitKind::Temperature {
                        // node, [11,] value
                        if k < toks.len() {
                            if parse_i32(&toks[k]).ok() == Some(11) {
                                k += 1;
                            }
                        }
                        if k >= toks.len() {
                            break;
                        }
                        let val = parse_f64(&toks[k])?;
                        k += 1;
                        if let Ok(id) = parse_i32(&name) {
                            model.init.push(InitCond {
                                node: id,
                                dof: 0,
                                value: val,
                                kind,
                            });
                        } else {
                            model.warnings.push(format!("__IC__|T|{name}|0|{val}"));
                        }
                    } else {
                        if k + 1 >= toks.len() {
                            break;
                        }
                        let dof = parse_i32(&toks[k])? as usize;
                        let val = parse_f64(&toks[k + 1])?;
                        k += 2;
                        if dof < 1 || dof > 6 {
                            continue;
                        }
                        if let Ok(id) = parse_i32(&name) {
                            model.init.push(InitCond {
                                node: id,
                                dof: dof - 1,
                                value: val,
                                kind,
                            });
                        } else {
                            let tag = match kind {
                                InitKind::Velocity => "V",
                                _ => "U",
                            };
                            model
                                .warnings
                                .push(format!("__IC__|{tag}|{name}|{}|{val}", dof - 1));
                        }
                    }
                }
            }
            "*EQUATION" => {
                if params.keys().any(|k| k.contains("REMOVE")) {
                    model.warn("*EQUATION, REMOVE wird ignoriert.");
                    let (_toks, ni) = collect_tokens(&lines, i + 1);
                    i = ni;
                    continue;
                }
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut k = 0;
                while k < toks.len() {
                    let nterms = parse_i32(&toks[k]).unwrap_or(0) as usize;
                    k += 1;
                    if nterms == 0 {
                        break;
                    }
                    let mut raw_terms: Vec<(String, usize, f64)> = Vec::new();
                    let mut rhs = 0.0;
                    for _ in 0..nterms {
                        if k + 2 >= toks.len() {
                            return err("*EQUATION: zu wenige Terme");
                        }
                        let node_tok = toks[k].to_ascii_uppercase();
                        let dof = parse_i32(&toks[k + 1])? as usize;
                        let coef = parse_f64(&toks[k + 2])?;
                        k += 3;
                        if dof >= 1 && dof <= 6 {
                            raw_terms.push((node_tok, dof - 1, coef));
                        }
                    }
                    if k < toks.len() && parse_i32(&toks[k]).is_err() {
                        if let Some(v) = parse_f64_inner(&toks[k]) {
                            rhs = v;
                            k += 1;
                        }
                    }
                    let named = raw_terms.iter().any(|(n, _, _)| parse_i32(n).is_err());
                    if named {
                        let payload: Vec<String> = raw_terms
                            .iter()
                            .map(|(n, d, c)| format!("{n}|{d}|{c}"))
                            .collect();
                        model
                            .warnings
                            .push(format!("__EQ__|{}|{}|{rhs}", raw_terms.len(), payload.join("|")));
                    } else {
                        let mut terms = Vec::new();
                        for (n, d, c) in raw_terms {
                            terms.push((parse_i32(&n)?, d, c));
                        }
                        model.equations.push(Equation { terms, rhs });
                    }
                }
            }
            "*SURFACE" => {
                let name = params
                    .get("NAME")
                    .cloned()
                    .ok_or_else(|| crate::error::FemError("*SURFACE ohne NAME=".into()))?;
                let ty = params
                    .get("TYPE")
                    .cloned()
                    .unwrap_or_else(|| "ELEMENT".into());
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut surf = Surface {
                    name: name.clone(),
                    nodes: Vec::new(),
                    faces: Vec::new(),
                };
                if ty == "NODE" {
                    for t in &toks {
                        if let Ok(id) = parse_i32(t) {
                            surf.nodes.push(id);
                        } else {
                            model
                                .warnings
                                .push(format!("__SN__|{name}|{}", t.to_ascii_uppercase()));
                        }
                    }
                } else {
                    let mut k = 0;
                    while k < toks.len() {
                        let label = toks[k].to_ascii_uppercase();
                        k += 1;
                        let face = if k < toks.len() {
                            let f = toks[k].to_ascii_uppercase();
                            if f.starts_with('S') {
                                k += 1;
                                if f == "SPOS" || f == "SNEG" {
                                    1
                                } else {
                                    f.trim_start_matches('S').parse::<i32>().unwrap_or(1)
                                }
                            } else {
                                1
                            }
                        } else {
                            1
                        };
                        if let Ok(id) = parse_i32(&label) {
                            surf.faces.push((id, face));
                        } else {
                            model
                                .warnings
                                .push(format!("__SF__|{name}|{label}|{face}"));
                        }
                    }
                }
                model.surfaces.insert(name.to_ascii_uppercase(), surf);
            }
            "*SURFACE INTERACTION" => {
                let name = params
                    .get("NAME")
                    .cloned()
                    .ok_or_else(|| crate::error::FemError("*SURFACE INTERACTION ohne NAME=".into()))?;
                let key = name.to_ascii_uppercase();
                model.interactions.insert(
                    key.clone(),
                    SurfaceInteraction {
                        kn: 1.0e7,
                        mu: 0.0,
                    },
                );
                current_interaction = Some(key);
                let (_toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
            }
            "*SURFACE BEHAVIOR" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let kn = toks.first().and_then(|t| parse_f64(t).ok()).unwrap_or(1.0e7);
                if let Some(name) = &current_interaction {
                    if let Some(it) = model.interactions.get_mut(name) {
                        it.kn = kn.abs().max(0.0);
                    }
                } else {
                    model.warn("*SURFACE BEHAVIOR ohne *SURFACE INTERACTION — ignoriert.");
                }
            }
            "*FRICTION" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mu = toks.first().and_then(|t| parse_f64(t).ok()).unwrap_or(0.0);
                if let Some(name) = &current_interaction {
                    if let Some(it) = model.interactions.get_mut(name) {
                        it.mu = mu.abs();
                    }
                }
            }
            "*CONTACT PAIR" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if toks.len() < 2 {
                    model.warn("*CONTACT PAIR braucht SLAVE, MASTER.");
                    continue;
                }
                let iname = params
                    .get("INTERACTION")
                    .cloned()
                    .unwrap_or_default()
                    .to_ascii_uppercase();
                let (kn, mu) = model
                    .interactions
                    .get(&iname)
                    .map(|it| (it.kn, it.mu))
                    .unwrap_or((1.0e7, 0.0));
                model.contact_pairs.push(ContactPair {
                    slave: toks[0].to_ascii_uppercase(),
                    master: toks[1].to_ascii_uppercase(),
                    kn,
                    mu,
                });
            }
            "*TIE" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let tol = params
                    .get("POSITION TOLERANCE")
                    .or_else(|| params.get("POSITIONTOLERANCE"))
                    .and_then(|s| parse_f64(s).ok())
                    .unwrap_or(0.0);
                if toks.len() >= 2 {
                    model.ties.push(Tie {
                        slave: toks[0].to_ascii_uppercase(),
                        master: toks[1].to_ascii_uppercase(),
                        position_tol: tol,
                    });
                }
            }
            "*RIGID BODY" => {
                i += 1;
                let nset = params
                    .get("NSET")
                    .cloned()
                    .ok_or_else(|| crate::error::FemError("*RIGID BODY ohne NSET=".into()))?;
                let refn = params
                    .get("REF NODE")
                    .or_else(|| params.get("REFNODE"))
                    .ok_or_else(|| crate::error::FemError("*RIGID BODY ohne REF NODE=".into()))?;
                let ref_node = parse_i32(refn)?;
                ensure_node(&mut model, ref_node);
                let mut rot_node = None;
                if let Some(rot) = params.get("ROT NODE").or_else(|| params.get("ROTNODE")) {
                    if let Ok(rid) = parse_i32(rot) {
                        ensure_node(&mut model, rid);
                        rot_node = Some(rid);
                    }
                }
                model.rigid_bodies.push(RigidBody {
                    nset: nset.to_ascii_uppercase(),
                    ref_node,
                    rot_node,
                });
            }
            "*COUPLING" => {
                let surf = params
                    .get("SURFACE")
                    .cloned()
                    .unwrap_or_default()
                    .to_ascii_uppercase();
                let refn = params
                    .get("REF NODE")
                    .or_else(|| params.get("REFNODE"))
                    .ok_or_else(|| crate::error::FemError("*COUPLING ohne REF NODE=".into()))?;
                let ref_node = parse_i32(refn)?;
                ensure_node(&mut model, ref_node);
                model.couplings.push(Coupling {
                    ref_node,
                    surface: surf,
                    kinematic: false,
                    dofs: vec![0, 1, 2],
                });
                current_coupling = Some(model.couplings.len() - 1);
                i += 1;
            }
            "*DISTRIBUTING" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if let Some(idx) = current_coupling {
                    model.couplings[idx].kinematic = false;
                    if toks.len() >= 2 {
                        let a = parse_i32(&toks[0]).unwrap_or(1);
                        let b = parse_i32(&toks[1]).unwrap_or(3);
                        model.couplings[idx].dofs = (a.max(1)..=b.max(a)).map(|d| (d - 1) as usize).collect();
                    }
                }
            }
            "*KINEMATIC" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if let Some(idx) = current_coupling {
                    model.couplings[idx].kinematic = true;
                    if toks.len() >= 2 {
                        let a = parse_i32(&toks[0]).unwrap_or(1);
                        let b = parse_i32(&toks[1]).unwrap_or(3);
                        model.couplings[idx].dofs = (a.max(1)..=b.max(a)).map(|d| (d - 1) as usize).collect();
                    }
                }
            }
            "*TRANSFORM" => {
                let nset = params
                    .get("NSET")
                    .cloned()
                    .ok_or_else(|| crate::error::FemError("*TRANSFORM ohne NSET=".into()))?;
                let cylindrical = params
                    .get("TYPE")
                    .map(|s| s.starts_with('C'))
                    .unwrap_or(false);
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if toks.len() < 6 {
                    return err("*TRANSFORM erwartet zwei Richtungsvektoren (6 Zahlen).");
                }
                let a = [
                    parse_f64(&toks[0])?,
                    parse_f64(&toks[1])?,
                    parse_f64(&toks[2])?,
                ];
                let b = [
                    parse_f64(&toks[3])?,
                    parse_f64(&toks[4])?,
                    parse_f64(&toks[5])?,
                ];
                let na = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
                let (e1, e2, e3, origin, axis) = if cylindrical || na < 1e-18 {
                    let axis = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                    let e3 = unit(axis);
                    let helper = if e3[2].abs() < 0.9 {
                        [0.0, 0.0, 1.0]
                    } else {
                        [1.0, 0.0, 0.0]
                    };
                    let mut e1 = cross(helper, e3);
                    let n1 = (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]).sqrt();
                    if n1 < 1e-18 {
                        e1 = [1.0, 0.0, 0.0];
                    } else {
                        e1 = [e1[0] / n1, e1[1] / n1, e1[2] / n1];
                    }
                    let e2 = cross(e3, e1);
                    (e1, e2, e3, a, axis)
                } else {
                    let e1 = unit(a);
                    let mut e3 = cross(e1, b);
                    let n3 = (e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2]).sqrt();
                    if n3 < 1e-18 {
                        let helper = if e1[2].abs() < 0.9 {
                            [0.0, 0.0, 1.0]
                        } else {
                            [1.0, 0.0, 0.0]
                        };
                        e3 = cross(e1, helper);
                        let n3b = (e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2]).sqrt().max(1e-30);
                        e3 = [e3[0] / n3b, e3[1] / n3b, e3[2] / n3b];
                        model.warn("*TRANSFORM: Achsen waren parallel — Hilfsachse verwendet.");
                    } else {
                        e3 = [e3[0] / n3, e3[1] / n3, e3[2] / n3];
                    }
                    let e2 = cross(e3, e1);
                    (e1, e2, e3, [0.0, 0.0, 0.0], e3)
                };
                let axes = [
                    [e1[0], e2[0], e3[0]],
                    [e1[1], e2[1], e3[1]],
                    [e1[2], e2[2], e3[2]],
                ];
                model.transforms.push(Transform {
                    nset: nset.to_ascii_uppercase(),
                    axes,
                    cylindrical: cylindrical || na < 1e-18,
                    origin,
                    axis,
                });
            }
            "*SPRING" => {
                let elset = params
                    .get("ELSET")
                    .cloned()
                    .unwrap_or_else(|| "EALL".into());
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                // optional dof line then stiffness, or just stiffness
                let k = if toks.len() >= 2 {
                    parse_f64(&toks[toks.len() - 1])?
                } else if toks.len() == 1 {
                    parse_f64(&toks[0])?
                } else {
                    return err("*SPRING ohne Steifigkeit");
                };
                model.elset_spring.insert(elset, k);
            }
            "*MASS" => {
                let elset = params
                    .get("ELSET")
                    .cloned()
                    .unwrap_or_else(|| "EALL".into());
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let m = if !toks.is_empty() {
                    parse_f64(&toks[0])?
                } else {
                    return err("*MASS ohne Wert");
                };
                model.elset_mass.insert(elset, m);
            }
            "*ROTARY INERTIA" | "*ROTARYI" => {
                let elset = params
                    .get("ELSET")
                    .cloned()
                    .unwrap_or_else(|| "EALL".into());
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut ijk = [0.0; 6];
                for (k, t) in toks.iter().take(6).enumerate() {
                    ijk[k] = parse_f64(t).unwrap_or(0.0);
                }
                model.elset_rotary.insert(elset, ijk);
            }
            "*DASHPOT" => {
                let elset = params
                    .get("ELSET")
                    .cloned()
                    .unwrap_or_else(|| "EALL".into());
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let c = if !toks.is_empty() {
                    parse_f64(&toks[toks.len() - 1])?
                } else {
                    return err("*DASHPOT ohne Dämpfung");
                };
                model.elset_dashpot.insert(elset, c);
            }
            "*GAP" => {
                let elset = params
                    .get("ELSET")
                    .cloned()
                    .unwrap_or_else(|| "EALL".into());
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let clearance = toks.get(0).and_then(|t| parse_f64(t).ok()).unwrap_or(0.0);
                let k = toks.get(1).and_then(|t| parse_f64(t).ok()).unwrap_or(1.0e8);
                let mut dir = [0.0; 3];
                if toks.len() >= 5 {
                    dir[0] = parse_f64(&toks[toks.len() - 3]).unwrap_or(0.0);
                    dir[1] = parse_f64(&toks[toks.len() - 2]).unwrap_or(0.0);
                    dir[2] = parse_f64(&toks[toks.len() - 1]).unwrap_or(0.0);
                }
                model.elset_gap.insert(
                    elset,
                    crate::model::GapSection { clearance, k, dir },
                );
            }
            "*HYPERELASTIC" | "*HYPERFOAM" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let name = current_material
                    .clone()
                    .unwrap_or_else(|| "MATERIAL-1".into());
                let m = model.materials.entry(name).or_default();
                let vals: Vec<f64> = toks.iter().filter_map(|t| parse_f64_inner(t)).collect();
                if kw == "*HYPERFOAM" {
                    m.hyper = HyperKind::Hyperfoam;
                    let n = params
                        .get("N")
                        .and_then(|s| s.parse::<u8>().ok())
                        .unwrap_or(1)
                        .clamp(1, 2);
                    m.h_n = n;
                    for (i, v) in vals.iter().copied().take(6).enumerate() {
                        m.h[i] = v;
                    }
                } else {
                    let typ = params
                        .keys()
                        .find(|k| k.contains("NEO") || k.contains("OGDEN") || k.contains("MOONEY"))
                        .cloned()
                        .unwrap_or_default();
                    if typ.contains("OGDEN") {
                        m.hyper = HyperKind::Ogden;
                        let n = params
                            .get("N")
                            .and_then(|s| s.parse::<u8>().ok())
                            .unwrap_or(1)
                            .clamp(1, 2);
                        m.h_n = n;
                        for (i, v) in vals.iter().copied().take(6).enumerate() {
                            m.h[i] = v;
                        }
                    } else {
                        m.hyper = HyperKind::NeoHooke;
                        m.h_n = 1;
                        m.h[0] = vals.first().copied().unwrap_or(1.0);
                        m.h[1] = vals.get(1).copied().unwrap_or(0.0);
                    }
                }
                m.set_equiv_from_hyper();
            }
            "*PLASTIC" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let name = current_material
                    .clone()
                    .unwrap_or_else(|| "MATERIAL-1".into());
                let mut pts = Vec::new();
                let mut k = 0;
                while k + 1 < toks.len() {
                    let sy = parse_f64(&toks[k])?;
                    let pe = parse_f64(&toks[k + 1])?;
                    k += 2;
                    pts.push((pe, sy));
                }
                pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                model.plastic.insert(name, pts);
            }
            "*EXPANSION" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if toks.is_empty() {
                    return err("*EXPANSION ohne Wert");
                }
                let alpha = parse_f64(&toks[0])?;
                let name = current_material
                    .clone()
                    .unwrap_or_else(|| "MATERIAL-1".into());
                model.materials.entry(name).or_default().alpha = alpha;
            }
            "*TEMPERATURE" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut k = 0;
                while k + 1 < toks.len() {
                    let name = toks[k].to_ascii_uppercase();
                    let t = parse_f64(&toks[k + 1])?;
                    k += 2;
                    model.warnings.push(format!("__TP__|{name}|{t}"));
                }
            }
            "*NODE FILE" | "*NODE OUTPUT" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                model.output_u = false;
                model.output_rf = false;
                for t in toks {
                    match t.to_ascii_uppercase().as_str() {
                        "U" => model.output_u = true,
                        "RF" => model.output_rf = true,
                        "NT" => model.output_nt = true,
                        _ => {}
                    }
                }
                if !model.output_u && !model.output_rf {
                    model.output_u = true;
                }
            }
            "*EL FILE" | "*ELEMENT OUTPUT" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                model.output_s = false;
                model.output_e = false;
                for t in toks {
                    match t.to_ascii_uppercase().as_str() {
                        "S" => model.output_s = true,
                        "E" => model.output_e = true,
                        _ => {}
                    }
                }
                if !model.output_s && !model.output_e {
                    model.output_s = true;
                }
            }
            "*NODE PRINT" | "*EL PRINT" => {
                let (_toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
            }
            "*INCLUDE" => {
                i += 1;
            }
            "*PREPRINT" => {
                i += 1;
            }
            "*CONTACT" => {
                let (_toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
            }
            "*END STEP" | "*END STEP " => {
                if step_open {
                    push_step(
                        &mut model,
                        step_bc_from,
                        step_cload_from,
                        step_dload_from,
                    );
                }
                step_open = false;
                i += 1;
            }
            "*CONTROLS" => {
                if let Some(m) = params.get("MAXITER") {
                    if let Ok(v) = parse_i32(m) {
                        model.max_newton = v.max(1) as usize;
                    }
                }
                if let Some(t) = params.get("RTOL").or_else(|| params.get("FINT")) {
                    if let Ok(v) = parse_f64(t) {
                        if v > 0.0 {
                            model.newton_tol = v;
                        }
                    }
                }
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if model.max_newton == 25 && !toks.is_empty() {
                    if let Ok(v) = parse_i32(&toks[0]) {
                        if v > 0 {
                            model.max_newton = v as usize;
                        }
                    }
                }
            }
            other => {
                if other.starts_with('*') {
                    match other {
                        "*COUPLED TEMPERATURE-DISPLACEMENT" | "*VISCO" | "*CREEP" => {
                            return err(format!(
                                "{other} nicht unterstützt in dieser Version."
                            ));
                        }
                        "*ORIENTATION" => {
                            model.warn(format!("{other} wird ignoriert."));
                            let (_toks, ni) = collect_tokens(&lines, i + 1);
                            i = ni;
                        }
                        _ => {
                            model.warn(format!("Unbekanntes Schlüsselwort {other} — ignoriert."));
                            let (_toks, ni) = collect_tokens(&lines, i + 1);
                            i = ni;
                        }
                    }
                } else {
                    i += 1;
                }
            }
        }
    }

    if model.node_ids.is_empty() {
        return err("Keine Knoten (*NODE) im Eingabedeck.");
    }
    if model.elements.is_empty() {
        return err("Keine Elemente (*ELEMENT) im Eingabedeck.");
    }
    if model.materials.is_empty() {
        model.warn("Kein *MATERIAL — E=210000, nu=0.3 angenommen.");
        model.materials.insert("STEEL".into(), Material::default());
    }
    if !saw_step {
        model.warn("Kein *STEP — linear-statischer Schritt angenommen.");
    }
    if step_open || model.steps.is_empty() {
        push_step(
            &mut model,
            step_bc_from,
            step_cload_from,
            step_dload_from,
        );
    }
    if let Some(last) = model.steps.last() {
        model.procedure = last.procedure.clone();
    }

    model.compact();
    expand_deferred(&mut model)?;
    if let Some(last) = model.steps.last_mut() {
        last.n_cload = last.n_cload.max(model.cloads.len());
        last.n_dload = last.n_dload.max(model.dloads.len());
        last.n_bc = last.n_bc.max(model.bcs.len());
    }
    Ok(model)
}

fn op_is_new(params: &HashMap<String, String>) -> bool {
    params
        .get("OP")
        .map(|v| v.eq_ignore_ascii_case("NEW"))
        .unwrap_or(false)
}

fn push_step(
    model: &mut Model,
    bc_from: usize,
    cload_from: usize,
    dload_from: usize,
) {
    model.steps.push(AnalysisStep {
        procedure: model.procedure.clone(),
        n_cload: model.cloads.len(),
        n_dload: model.dloads.len(),
        n_bc: model.bcs.len(),
        cload_from,
        dload_from,
        bc_from,
    });
}

fn parse_beam_section(sectyp: &str, toks: &[String]) -> Result<BeamSection> {
    let ty = sectyp.to_ascii_uppercase();
    let mut nums = Vec::new();
    for t in toks {
        nums.push(parse_f64(t)?);
    }
    let n_dim = match ty.as_str() {
        "RECT" | "RECTANGULAR" => 2,
        "CIRC" | "CIRCULAR" => 1,
        "PIPE" => 2,
        "BOX" => {
            if nums.len() >= 9 {
                6
            } else if nums.len() >= 7 {
                4
            } else {
                2
            }
        }
        "GENERAL" | "ARBITRARY" => 5,
        other => {
            return err(format!(
                "BEAM SECTION={other} nicht unterstützt. RECT, CIRC, PIPE, BOX, GENERAL."
            ));
        }
    };
    if nums.len() < n_dim {
        return err(format!(
            "*BEAM SECTION {ty} erwartet mindestens {n_dim} Werte"
        ));
    }
    let mut n1 = [0.0, 0.0, -1.0];
    if nums.len() >= n_dim + 3 {
        n1 = [nums[n_dim], nums[n_dim + 1], nums[n_dim + 2]];
        if (n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]).sqrt() < 1e-18 {
            n1 = [0.0, 0.0, -1.0];
        }
    }
    Ok(match ty.as_str() {
        "RECT" | "RECTANGULAR" => BeamSection::rect(nums[0], nums[1], n1),
        "CIRC" | "CIRCULAR" => BeamSection::circ(nums[0], n1),
        "PIPE" => BeamSection::pipe(nums[0], nums[1], n1),
        "BOX" => {
            let t = if nums.len() >= 3 { nums[2] } else { nums[0].min(nums[1]) * 0.1 };
            let t2 = if nums.len() >= 4 { nums[3] } else { t };
            let t3 = if nums.len() >= 5 { nums[4] } else { t };
            let t4 = if nums.len() >= 6 { nums[5] } else { t };
            BeamSection::box_sec(nums[0], nums[1], t, t2, t3, t4, n1)
        }
        _ => BeamSection::general(nums[0], nums[1], nums[2], nums[3], nums[4], n1),
    })
}

fn expand_deferred(model: &mut Model) -> Result<()> {
    let warnings = std::mem::take(&mut model.warnings);
    let mut nested_e: Vec<(String, String)> = Vec::new();
    let mut nested_n: Vec<(String, String)> = Vec::new();
    let mut rest = Vec::new();
    for w in warnings {
        if let Some(r) = w.strip_prefix("__ES__|") {
            let p: Vec<&str> = r.split('|').collect();
            if p.len() >= 2 {
                nested_e.push((p[0].to_string(), p[1].to_string()));
            }
        } else if let Some(r) = w.strip_prefix("__NS__|") {
            let p: Vec<&str> = r.split('|').collect();
            if p.len() >= 2 {
                nested_n.push((p[0].to_string(), p[1].to_string()));
            }
        } else {
            rest.push(w);
        }
    }
    expand_nested_sets(&mut model.elsets, &nested_e, "Elementset")?;
    expand_nested_sets(&mut model.nsets, &nested_n, "Knotenset")?;
    let mut keep_warn = Vec::new();
    let warnings = rest;
    for w in warnings {
        if let Some(rest) = w.strip_prefix("__BC__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 4 {
                let name = p[0];
                let dof: usize = p[1].parse().unwrap_or(0);
                let val: f64 = p[3].parse().unwrap_or(0.0);
                let nodes = model.expand_nset(name)?;
                for node in nodes {
                    model.bcs.push(Boundary {
                        node,
                        dof,
                        value: val,
                    });
                }
            }
        } else if let Some(rest) = w.strip_prefix("__CL__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 3 {
                let name = p[0];
                let dof: usize = p[1].parse().unwrap_or(0);
                let mag: f64 = p[2].parse().unwrap_or(0.0);
                let nodes = model.expand_nset(name)?;
                for node in nodes {
                    model.cloads.push(Cload {
                        node,
                        dof,
                        mag,
                        amplitude: if p.len() >= 4 {
                            p[3].to_string()
                        } else {
                            String::new()
                        },
                    });
                }
            }
        } else if let Some(rest) = w.strip_prefix("__DL__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 3 {
                let name = p[0];
                let face: i32 = p[1].parse().unwrap_or(0);
                let mag: f64 = p[2].parse().unwrap_or(0.0);
                let elems = model.expand_elset(name)?;
                for elem in elems {
                    model.dloads.push(Dload::Pressure { elem, face, mag });
                }
            }
        } else if let Some(rest) = w.strip_prefix("__BG__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 5 {
                let name = p[0];
                let dx: f64 = p[1].parse().unwrap_or(0.0);
                let dy: f64 = p[2].parse().unwrap_or(0.0);
                let dz: f64 = p[3].parse().unwrap_or(0.0);
                let mag: f64 = p[4].parse().unwrap_or(0.0);
                let elems = model.expand_elset(name)?;
                for elem in elems {
                    model.dloads.push(Dload::BeamGlobal {
                        elem,
                        dir: [dx, dy, dz],
                        mag,
                    });
                }
            }
        } else if let Some(rest) = w.strip_prefix("__CFUG__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 8 {
                let name = p[0];
                let omega2: f64 = p[1].parse().unwrap_or(0.0);
                let p1 = [
                    p[2].parse().unwrap_or(0.0),
                    p[3].parse().unwrap_or(0.0),
                    p[4].parse().unwrap_or(0.0),
                ];
                let p2 = [
                    p[5].parse().unwrap_or(0.0),
                    p[6].parse().unwrap_or(0.0),
                    p[7].parse().unwrap_or(0.0),
                ];
                let elems = model.expand_elset(name).unwrap_or_default();
                model.dloads.push(Dload::Centrif {
                    omega2,
                    p1,
                    p2,
                    elems,
                });
            }
        } else if let Some(rest) = w.strip_prefix("__EQ__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 5 {
                let nterms: usize = p[0].parse().unwrap_or(0);
                let rhs: f64 = p.last().and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let mut raw: Vec<(String, usize, f64)> = Vec::new();
                let mut j = 1;
                for _ in 0..nterms {
                    if j + 2 >= p.len() {
                        break;
                    }
                    let name = p[j].to_string();
                    let dof: usize = p[j + 1].parse().unwrap_or(0);
                    let coef: f64 = p[j + 2].parse().unwrap_or(0.0);
                    raw.push((name, dof, coef));
                    j += 3;
                }
                let named_idx = raw.iter().position(|(n, _, _)| parse_i32(n).is_err());
                if let Some(si) = named_idx {
                    if let Ok(nodes) = model.expand_nset(&raw[si].0) {
                        for node in nodes {
                            let mut terms = Vec::new();
                            let mut ok = true;
                            for (idx, (n, d, c)) in raw.iter().enumerate() {
                                let nid = if idx == si {
                                    node
                                } else if let Ok(id) = parse_i32(n) {
                                    id
                                } else if let Ok(more) = model.expand_nset(n) {
                                    if more.len() == 1 {
                                        more[0]
                                    } else {
                                        ok = false;
                                        break;
                                    }
                                } else {
                                    ok = false;
                                    break;
                                };
                                terms.push((nid, *d, *c));
                            }
                            if ok {
                                model.equations.push(Equation { terms, rhs });
                            }
                        }
                    }
                }
            }
        } else if let Some(rest) = w.strip_prefix("__DLN__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 3 {
                let name = p[0];
                let face: i32 = p[1].parse().unwrap_or(1);
                let nid: i32 = p[2].parse().unwrap_or(0);
                let mag = model
                    .bcs
                    .iter()
                    .find(|b| b.node == nid)
                    .map(|b| b.value)
                    .unwrap_or(1.0);
                let elems = if let Ok(id) = parse_i32(name) {
                    vec![id]
                } else {
                    model.expand_elset(name).unwrap_or_default()
                };
                for elem in elems {
                    model.dloads.push(Dload::Pressure { elem, face, mag });
                }
            }
        } else if let Some(rest) = w.strip_prefix("__SN__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 2 {
                let sname = p[0].to_ascii_uppercase();
                let set = p[1];
                if let Ok(nodes) = model.expand_nset(set) {
                    if let Some(s) = model.surfaces.get_mut(&sname) {
                        s.nodes.extend(nodes);
                    }
                }
            }
        } else if let Some(rest) = w.strip_prefix("__SF__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 3 {
                let sname = p[0].to_ascii_uppercase();
                let set = p[1];
                let face: i32 = p[2].parse().unwrap_or(1);
                if let Ok(elems) = model.expand_elset(set) {
                    if let Some(s) = model.surfaces.get_mut(&sname) {
                        for e in elems {
                            s.faces.push((e, face));
                        }
                    }
                }
            }
        } else if let Some(rest) = w.strip_prefix("__TP__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 2 {
                let name = p[0];
                let t: f64 = p[1].parse().unwrap_or(0.0);
                if let Ok(nodes) = model.expand_nset(name) {
                    for n in nodes {
                        model.temperatures.insert(n, t);
                    }
                }
            }
        } else if let Some(rest) = w.strip_prefix("__TB__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 2 {
                let name = p[0];
                let val: f64 = p[1].parse().unwrap_or(0.0);
                if let Ok(nodes) = model.expand_nset(name) {
                    for node in nodes {
                        model.thermal_bcs.push(ThermalBc { node, value: val });
                    }
                }
            }
        } else if let Some(rest) = w.strip_prefix("__CF__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 2 {
                let name = p[0];
                let mag: f64 = p[1].parse().unwrap_or(0.0);
                if let Ok(nodes) = model.expand_nset(name) {
                    for node in nodes {
                        model.cfluxes.push(Cflux { node, mag });
                    }
                }
            }
        } else if let Some(rest) = w.strip_prefix("__DF__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 3 {
                let name = p[0];
                let tag: i32 = p[1].parse().unwrap_or(0);
                let mag: f64 = p[2].parse().unwrap_or(0.0);
                let kind = if tag == 0 {
                    FluxKind::Body
                } else {
                    FluxKind::Face(tag)
                };
                if let Ok(elems) = model.expand_elset(name) {
                    for elem in elems {
                        model.dfluxes.push(Dflux { elem, kind, mag });
                    }
                }
            }
        } else if let Some(rest) = w.strip_prefix("__FM__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 4 {
                let name = p[0];
                let face: i32 = p[1].parse().unwrap_or(1);
                let tinf: f64 = p[2].parse().unwrap_or(0.0);
                let h: f64 = p[3].parse().unwrap_or(0.0);
                if let Ok(elems) = model.expand_elset(name) {
                    for elem in elems {
                        model.films.push(Film {
                            elem,
                            face,
                            t_inf: tinf,
                            h,
                        });
                    }
                }
            }
        } else if let Some(rest) = w.strip_prefix("__IC__|") {
            let p: Vec<&str> = rest.split('|').collect();
            if p.len() >= 4 {
                let kind = match p[0] {
                    "V" => InitKind::Velocity,
                    "U" => InitKind::Displacement,
                    _ => InitKind::Temperature,
                };
                let name = p[1];
                let dof: usize = p[2].parse().unwrap_or(0);
                let val: f64 = p[3].parse().unwrap_or(0.0);
                if let Ok(nodes) = model.expand_nset(name) {
                    for node in nodes {
                        model.init.push(InitCond {
                            node,
                            dof,
                            value: val,
                            kind,
                        });
                    }
                }
            }
        } else {
            keep_warn.push(w);
        }
    }
    model.warnings = keep_warn;
    Ok(())
}

fn expand_nested_sets(
    sets: &mut HashMap<String, Vec<i32>>,
    refs: &[(String, String)],
    kind: &str,
) -> Result<()> {
    if refs.is_empty() {
        return Ok(());
    }
    let mut pending = refs.to_vec();
    for _ in 0..16 {
        if pending.is_empty() {
            return Ok(());
        }
        let before = pending.len();
        let mut next = Vec::new();
        for (parent, child) in pending {
            if let Ok(id) = child.parse::<i32>() {
                sets.entry(parent).or_default().push(id);
                continue;
            }
            if let Some(ids) = sets.get(&child).cloned() {
                sets.entry(parent).or_default().extend(ids);
            } else {
                next.push((parent, child));
            }
        }
        pending = next;
        if pending.len() == before {
            break;
        }
    }
    if !pending.is_empty() {
        return err(format!(
            "Unbekanntes {kind} {} (in {})",
            pending[0].1, pending[0].0
        ));
    }
    Ok(())
}

fn unit(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-30);
    [v[0] / n, v[1] / n, v[2] / n]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
