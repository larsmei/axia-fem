use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::error::{err, Result};
use crate::model::{
    BeamSection, Boundary, Cload, Coupling, Dload, ElemKind, Element, Equation, Material, Model,
    RigidBody, Surface, Tie, Transform,
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
    (kw, params)
}

fn tokenize_data(line: &str) -> Vec<String> {
    line.split(|c: char| c == ',' || c.is_whitespace())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect()
}

fn parse_f64(s: &str) -> Result<f64> {
    let t = s.trim().trim_end_matches('.');
    t.parse::<f64>()
        .or_else(|_| s.trim().parse::<f64>())
        .map_err(|_| crate::error::FemError(format!("Keine Zahl: {s}")))
}

fn parse_i32(s: &str) -> Result<i32> {
    let t = s.trim();
    if let Ok(v) = t.parse::<i32>() {
        return Ok(v);
    }
    let f = parse_f64(t)?;
    Ok(f.round() as i32)
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
    let mut saw_step = false;

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
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                if toks.is_empty() {
                    continue;
                }
                let stride = if toks.len() % 4 == 0 {
                    4
                } else if toks.len() % 3 == 0 {
                    3
                } else if toks.len() >= 4 {
                    4
                } else {
                    3
                };
                let mut k = 0;
                let mut added = Vec::new();
                while k + stride - 1 < toks.len() {
                    let id = parse_i32(&toks[k])?;
                    let x = parse_f64(&toks[k + 1])?;
                    let y = parse_f64(&toks[k + 2])?;
                    let z = if stride == 4 {
                        parse_f64(&toks[k + 3])?
                    } else {
                        0.0
                    };
                    k += stride;
                    model.node_ids.push(id);
                    model.coords.push([x, y, z]);
                    added.push(id);
                }
                if k != toks.len() {
                    model.warn(format!(
                        "*NODE: {} Zahlen übrig, ignoriert.",
                        toks.len() - k
                    ));
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
                        "Nicht unterstützter Elementtyp {typ}. Unterstützt: C3D8/C3D8I/C3D8R, C3D20, C3D4, C3D10, C3D6, CPS*, CPE*, S3/S4/S6/S8, B31/B32, T3D2/T3D3, SPRINGA."
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
                        ids.push(parse_i32(t)?);
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
                        ids.push(parse_i32(t)?);
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
                let e = parse_f64(&toks[0])?;
                let nu = parse_f64(&toks[1])?;
                let name = current_material
                    .clone()
                    .unwrap_or_else(|| "MATERIAL-1".into());
                let m = model.materials.entry(name).or_default();
                m.e = e;
                m.nu = nu;
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
            "*SOLID SECTION" | "*SHELL SECTION" => {
                let elset = params.get("ELSET").cloned().unwrap_or_else(|| "EALL".into());
                if let Some(mat) = params.get("MATERIAL") {
                    model.elset_material.insert(elset.clone(), mat.clone());
                }
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let t = if !toks.is_empty() {
                    parse_f64(&toks[0]).unwrap_or(1.0)
                } else {
                    1.0
                };
                model.elset_thickness.insert(elset, t);
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
                    let val = if k < toks.len() {
                        // Could be next set name or a value. If it parses as float and
                        // the token isn't an integer set name used later... values are
                        // typically 0 or a float. If the next-next looks like a dof, this
                        // is a new card.
                        // Cards: name, d1, [d2, [val]]
                        // Heuristic: if token has a dot/e or is 0, treat as value; if
                        // remaining tokens after taking it wouldn't form a new card, take it.
                        let peek = &toks[k];
                        if peek.contains('.') || peek.contains('e') || peek.contains('E') {
                            let v = parse_f64(peek)?;
                            k += 1;
                            v
                        } else if parse_i32(peek).is_ok() && k + 1 < toks.len() {
                            // Could be value 0 or start of next (node id).
                            // If following token is a dof 1-6, this is a new card.
                            let maybe_dof = parse_i32(&toks[k + 1]).ok();
                            if maybe_dof == Some(1)
                                || maybe_dof == Some(2)
                                || maybe_dof == Some(3)
                                || maybe_dof == Some(4)
                                || maybe_dof == Some(5)
                                || maybe_dof == Some(6)
                            {
                                0.0
                            } else {
                                let v = parse_f64(peek).unwrap_or(0.0);
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
                                // defer: keep a placeholder by recording the name in a special
                                // node id list via a side channel — we re-expand after compact.
                                Vec::new()
                            })
                            .collect()
                    };
                    // Store deferred BCs as node=-999999 and stash name in a warning-free list.
                    // We'll expand after all sets exist: keep raw cards.
                    if nodes.is_empty() && parse_i32(&name).is_err() {
                        // record as deferred using a sentinel with name in heading? Use a
                        // dedicated deferred buffer.
                        // We'll push a Boundary per node after a second pass.
                        // Encode the set name by pushing a fake node and keeping a parallel list.
                    }
                    let first = d1.max(1);
                    let last = d2.max(first);
                    if parse_i32(&name).is_ok() {
                        let id = parse_i32(&name)?;
                        for d in first..=last {
                            if d >= 1 && d <= 6 {
                                model.bcs.push(Boundary {
                                    node: id,
                                    dof: d - 1,
                                    value: val,
                                });
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
                            }
                        }
                    }
                }
            }
            "*CLOAD" => {
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
                        });
                    } else {
                        model.warnings.push(format!("__CL__|{name}|{}|{mag}", dof - 1));
                    }
                }
            }
            "*DLOAD" => {
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
                    } else if typ.starts_with('P') {
                        if k >= toks.len() {
                            return err("*DLOAD P ohne Betrag");
                        }
                        let mag = parse_f64(&toks[k])?;
                        k += 1;
                        let face = if typ == "P" {
                            0
                        } else {
                            typ.trim_start_matches('P').parse::<i32>().unwrap_or(0)
                        };
                        if let Ok(id) = parse_i32(&name) {
                            model.dloads.push(Dload::Pressure {
                                elem: id,
                                face,
                                mag,
                            });
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
                saw_step = true;
                if params.contains_key("NLGEOM") {
                    model.procedure = crate::model::Procedure::Static {
                        nlgeom: true,
                        increments: 1,
                    };
                    model.warn("NLGEOM: geometrisch nichtlineare Statik (falls unterstützt).");
                }
                i += 1;
            }
            "*STATIC" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let inc = if !toks.is_empty() {
                    parse_f64(&toks[0]).ok().map(|v| v.max(1.0) as usize).unwrap_or(1)
                } else {
                    1
                };
                let nlgeom = matches!(
                    model.procedure,
                    crate::model::Procedure::Static { nlgeom: true, .. }
                );
                model.procedure = crate::model::Procedure::Static {
                    nlgeom,
                    increments: inc.max(1),
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
            "*EQUATION" => {
                let (toks, ni) = collect_tokens(&lines, i + 1);
                i = ni;
                let mut k = 0;
                while k < toks.len() {
                    let nterms = parse_i32(&toks[k]).unwrap_or(0) as usize;
                    k += 1;
                    if nterms == 0 {
                        break;
                    }
                    let mut terms = Vec::new();
                    let mut rhs = 0.0;
                    for _ in 0..nterms {
                        if k + 2 >= toks.len() {
                            return err("*EQUATION: zu wenige Terme");
                        }
                        let node = parse_i32(&toks[k])?;
                        let dof = parse_i32(&toks[k + 1])? as usize;
                        let coef = parse_f64(&toks[k + 2])?;
                        k += 3;
                        if dof >= 1 && dof <= 6 {
                            terms.push((node, dof - 1, coef));
                        }
                    }
                    // optional constant
                    if k < toks.len() && parse_i32(&toks[k]).is_err() {
                        if let Ok(v) = parse_f64(&toks[k]) {
                            rhs = v;
                            k += 1;
                        }
                    }
                    model.equations.push(Equation { terms, rhs });
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
                model.rigid_bodies.push(RigidBody {
                    nset: nset.to_ascii_uppercase(),
                    ref_node,
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
                let e1 = unit(a);
                let mut e3 = cross(e1, b);
                let n3 = (e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2]).sqrt();
                if n3 < 1e-18 {
                    return err("*TRANSFORM: Vektoren sind parallel.");
                }
                e3[0] /= n3;
                e3[1] /= n3;
                e3[2] /= n3;
                let e2 = cross(e3, e1);
                // columns = local axes in global
                let axes = [
                    [e1[0], e2[0], e3[0]],
                    [e1[1], e2[1], e3[1]],
                    [e1[2], e2[2], e3[2]],
                ];
                model.transforms.push(Transform {
                    nset: nset.to_ascii_uppercase(),
                    axes,
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
            "*PREPRINT" | "*END STEP" | "*END STEP " => {
                i += 1;
            }
            other => {
                if other.starts_with('*') {
                    match other {
                        "*DYNAMIC" | "*HEAT TRANSFER"
                        | "*COUPLED TEMPERATURE-DISPLACEMENT" | "*VISCO" | "*CREEP" => {
                            return err(format!(
                                "{other} nicht unterstützt in dieser Version."
                            ));
                        }
                        "*PLASTIC" | "*CONTACT" | "*AMPLITUDE" | "*INITIAL CONDITIONS"
                        | "*ORIENTATION" => {
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

    model.compact();
    expand_deferred(&mut model)?;
    Ok(model)
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
        "GENERAL" | "ARBITRARY" => 5,
        other => {
            return err(format!(
                "BEAM SECTION={other} nicht unterstützt. RECT, CIRC, PIPE, GENERAL."
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
        _ => BeamSection::general(nums[0], nums[1], nums[2], nums[3], nums[4], n1),
    })
}

fn expand_deferred(model: &mut Model) -> Result<()> {
    let mut keep_warn = Vec::new();
    let warnings = std::mem::take(&mut model.warnings);
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
                    model.cloads.push(Cload { node, dof, mag });
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
        } else {
            keep_warn.push(w);
        }
    }
    model.warnings = keep_warn;
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
