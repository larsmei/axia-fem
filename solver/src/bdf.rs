//! MYSTRAN / Nastran bulk-data reader.
//!
//! Free-field and 8/16-character fixed fields, `$` comments, `INCLUDE`,
//! continuations, `SOL 1/101` statics and `SOL 3/103` modes. The result is an
//! Axia [`Model`](crate::model::Model). Shell `PLOAD2`/`PLOAD4` follow the
//! Nastran sign: positive pressure opposes the element normal.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use crate::error::{err, Result};
use crate::mystran_manifest::{self, CardStatus};
use crate::model::{
    AnalysisStep, BeamSection, Boundary, Cload, Dload, ElemKind, Element, Material, Model,
    Procedure, RigidBody,
};

pub fn is_mystran_deck(text: &str) -> bool {
    let mut bulk = false;
    let mut ccx = false;
    for line in text.lines() {
        let t = strip_dollar(line).trim();
        if t.is_empty() {
            continue;
        }
        let u = t.to_ascii_uppercase();
        if u.starts_with("BEGIN BULK") || u.starts_with("ENDDATA") || u == "CEND" || u.starts_with("CEND ")
        {
            bulk = true;
        }
        if u.starts_with('*') && !u.starts_with("**") {
            ccx = true;
        }
    }
    if bulk {
        return true;
    }
    if ccx {
        return false;
    }
    for line in text.lines() {
        let t = strip_dollar(line).trim();
        if t.is_empty() {
            continue;
        }
        let key = t
            .split([',', ' ', '\t'])
            .next()
            .unwrap_or("")
            .trim_end_matches('*')
            .to_ascii_uppercase();
        return matches!(
            key.as_str(),
            "GRID"
                | "GRDSET"
                | "CQUAD4"
                | "CTRIA3"
                | "CROD"
                | "CBAR"
                | "CBEAM"
                | "CONROD"
                | "CTETRA"
                | "CHEXA"
                | "CPENTA"
                | "MAT1"
                | "PSHELL"
                | "PSOLID"
                | "PROD"
                | "PBAR"
                | "PBARL"
                | "SPC"
                | "SPC1"
                | "FORCE"
                | "MOMENT"
                | "PLOAD2"
                | "PLOAD4"
                | "GRAV"
        );
    }
    false
}

/// Bulk names the parser accepts. Kept in lockstep with the manifest.
pub const RECOGNIZED_BULK: &[&str] = &[
    "PARAM", "DEBUG", "EIGRL", "GRDSET", "GRID", "CORD2R", "MAT1", "PSHELL", "PSOLID", "PROD",
    "PBAR", "PBARL", "CROD", "CONROD", "CBAR", "CBEAM", "CQUAD4", "CTRIA3", "CTETRA", "CHEXA",
    "CPENTA", "CELAS2", "CONM2", "RBE2", "FORCE", "MOMENT", "PLOAD2", "PLOAD4", "GRAV", "LOAD",
    "SPC", "SPC1", "SPCADD",
];

pub fn parse_with_base(text: &str, base: Option<&Path>) -> Result<Model> {
    let text = expand_includes(text, base, 0)?;
    let (exec, case_lines, bulk_lines) = split_sections(&text);
    let mut sol = 101i32;
    let mut id_title = String::new();
    for line in &exec {
        let (k, v) = split_kv(line);
        match k.as_str() {
            "SOL" => {
                sol = v
                    .split_whitespace()
                    .next()
                    .and_then(|s| parse_i32(s).ok())
                    .unwrap_or(101);
            }
            "ID" => {
                if id_title.is_empty() {
                    id_title = v;
                }
            }
            _ => {}
        }
    }
    let cases = parse_cases(&case_lines);
    let cards = assemble_cards(&bulk_lines)?;
    build_model(sol, id_title, &cases, &cards)
}

fn expand_includes(text: &str, base: Option<&Path>, depth: usize) -> Result<String> {
    if depth > 16 {
        return err("INCLUDE: Verschachtelung zu tief.");
    }
    let mut out = String::new();
    for line in text.lines() {
        let raw = strip_dollar(line).trim();
        let key = raw
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_ascii_uppercase();
        if key == "INCLUDE" {
            let spec = raw.split_once(char::is_whitespace).map(|(_, r)| r.trim()).unwrap_or("");
            let file = spec.trim_matches(|c| c == '\'' || c == '"').trim();
            if file.is_empty() {
                return err("INCLUDE ohne Dateiname.");
            }
            let path = match base {
                Some(b) => b.join(file),
                None => PathBuf::from(file),
            };
            #[cfg(target_arch = "wasm32")]
            {
                let _ = path;
                return err("INCLUDE ist im Browser nicht verfügbar.");
            }
            #[cfg(not(target_arch = "wasm32"))]
            {
                let nested = std::fs::read_to_string(&path).map_err(|e| {
                    crate::error::FemError(format!("INCLUDE kann '{}' nicht lesen: {e}", path.display()))
                })?;
                let nested_base = path.parent();
                out.push_str(&expand_includes(&nested, nested_base, depth + 1)?);
                if !out.ends_with('\n') {
                    out.push('\n');
                }
                continue;
            }
        }
        out.push_str(line);
        out.push('\n');
    }
    Ok(out)
}

fn split_sections(text: &str) -> (Vec<String>, Vec<String>, Vec<String>) {
    let mut exec = Vec::new();
    let mut case = Vec::new();
    let mut bulk = Vec::new();
    let mut phase = 0u8; // 0 exec, 1 case, 2 bulk
    let mut saw_begin = false;
    for line in text.lines() {
        let t = strip_dollar(line).trim();
        if t.is_empty() {
            continue;
        }
        let u = t.to_ascii_uppercase();
        if u.starts_with("BEGIN BULK") {
            phase = 2;
            saw_begin = true;
            continue;
        }
        if u.starts_with("ENDDATA") {
            break;
        }
        if phase == 0 && (u == "CEND" || u.starts_with("CEND ")) {
            phase = 1;
            continue;
        }
        match phase {
            0 => exec.push(t.to_string()),
            1 => case.push(line_keep_value(line)),
            _ => bulk.push(strip_dollar(line).trim_end().to_string()),
        }
    }
    if !saw_begin && bulk.is_empty() {
        // Bulk-only deck: everything that was classified as exec/case and is a card.
        let mut rest = exec;
        rest.extend(case.iter().map(|s| s.trim().to_string()));
        return (Vec::new(), Vec::new(), rest);
    }
    (exec, case, bulk)
}

/// Case-control value keeps the original spelling (titles).
fn line_keep_value(line: &str) -> String {
    strip_dollar(line).trim().to_string()
}

fn split_kv(line: &str) -> (String, String) {
    let t = line.trim();
    if let Some((a, b)) = t.split_once('=') {
        let key = a
            .split('(')
            .next()
            .unwrap_or(a)
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_ascii_uppercase();
        return (key, b.trim().to_string());
    }
    let mut it = t.split_whitespace();
    let key = it.next().unwrap_or("").to_ascii_uppercase();
    let val = it.collect::<Vec<_>>().join(" ");
    (key, val)
}

#[derive(Clone)]
struct CaseCtrl {
    title: String,
    subtitle: String,
    label: String,
    spc: Option<i32>,
    load: Option<i32>,
    method: Option<i32>,
}

fn parse_cases(lines: &[String]) -> Vec<CaseCtrl> {
    let mut global = CaseCtrl {
        title: String::new(),
        subtitle: String::new(),
        label: String::new(),
        spc: None,
        load: None,
        method: None,
    };
    let mut cases: Vec<CaseCtrl> = Vec::new();
    let mut cur: Option<CaseCtrl> = None;
    for line in lines {
        let (k, v) = split_kv(line);
        if k == "SUBCASE" {
            if let Some(c) = cur.take() {
                cases.push(c);
            }
            cur = Some(global.clone());
            continue;
        }
        let dest = cur.as_mut().unwrap_or(&mut global);
        match k.as_str() {
            "TITLE" => dest.title = v,
            "SUBTITLE" => dest.subtitle = v,
            "LABEL" => dest.label = v,
            "SPC" => dest.spc = v.split_whitespace().next().and_then(|s| parse_i32(s).ok()),
            "LOAD" => dest.load = v.split_whitespace().next().and_then(|s| parse_i32(s).ok()),
            "METHOD" => dest.method = v.split_whitespace().next().and_then(|s| parse_i32(s).ok()),
            // ECHO=SORT/UNSORT/NONE/BOTH changes the punch, not the solution.
            "ECHO" | "DISPLACEMENT" | "SPCFORCES" | "STRESS" | "FORCE" | "ELFORCE" | "OLOAD"
            | "STRAIN" | "MAXLINES" => {}
            _ => {}
        }
    }
    if let Some(c) = cur {
        cases.push(c);
    }
    if cases.is_empty() {
        cases.push(global);
    }
    cases
}

fn assemble_cards(lines: &[String]) -> Result<Vec<Vec<String>>> {
    let mut cards: Vec<Vec<String>> = Vec::new();
    let mut large = false;
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        if is_continuation(line) {
            if cards.is_empty() {
                return err(format!("Fortsetzungszeile ohne Karte: {line}"));
            }
            let extra = continuation_fields(line, large)?;
            cards.last_mut().unwrap().extend(extra);
            continue;
        }
        let (fields, is_large) = split_card(line)?;
        large = is_large;
        if fields.first().map(|s| s.is_empty()).unwrap_or(true) {
            continue;
        }
        cards.push(fields);
    }
    Ok(cards)
}

fn is_continuation(line: &str) -> bool {
    match line.as_bytes().first().copied() {
        Some(b'+') | Some(b'*') | Some(b',') | Some(b' ') => true,
        _ => false,
    }
}

fn split_card(line: &str) -> Result<(Vec<String>, bool)> {
    if line.contains(',') {
        let mut f: Vec<String> = line.split(',').map(|s| s.trim().to_string()).collect();
        if f.first().map(|s| s.is_empty()).unwrap_or(false) {
            return err(format!("Leerer Kartenname: {line}"));
        }
        f[0] = f[0].to_ascii_uppercase();
        while f.last().map(|s| s.is_empty()).unwrap_or(false) {
            f.pop();
        }
        return Ok((f, false));
    }
    let mut raw = line.to_string();
    if raw.len() < 80 {
        raw.push_str(&" ".repeat(80 - raw.len()));
    }
    let name_fld = raw[..8].trim().to_string();
    let large = name_fld.contains('*');
    let name = name_fld.trim_end_matches('*').trim().to_ascii_uppercase();
    let mut fields = vec![name];
    if large {
        for i in 0..4 {
            let a = 8 + i * 16;
            let b = (a + 16).min(raw.len());
            fields.push(raw.get(a..b).unwrap_or("").trim().to_string());
        }
    } else {
        for i in 1..9 {
            let a = i * 8;
            let b = a + 8;
            fields.push(raw.get(a..b).unwrap_or("").trim().to_string());
        }
    }
    Ok((fields, large))
}

fn continuation_fields(line: &str, large: bool) -> Result<Vec<String>> {
    if line.contains(',') || line.trim_start().starts_with(',') {
        let parts: Vec<String> = line.split(',').map(|s| s.trim().to_string()).collect();
        // Leading comma or "+ID" occupies the first token.
        let rest = if parts.first().map(|s| s.is_empty() || s.starts_with('+') || s.starts_with('*')).unwrap_or(false)
        {
            parts.into_iter().skip(1).collect::<Vec<_>>()
        } else {
            parts
        };
        return Ok(rest);
    }
    let mut raw = line.to_string();
    if raw.len() < 80 {
        raw.push_str(&" ".repeat(80 - raw.len()));
    }
    let mut fields = Vec::new();
    if large || raw[..8].contains('*') {
        for i in 0..4 {
            let a = 8 + i * 16;
            let b = (a + 16).min(raw.len());
            fields.push(raw.get(a..b).unwrap_or("").trim().to_string());
        }
    } else {
        for i in 1..9 {
            let a = i * 8;
            fields.push(raw.get(a..a + 8).unwrap_or("").trim().to_string());
        }
    }
    Ok(fields)
}

#[derive(Clone)]
struct MatRec {
    e: f64,
    nu: f64,
    rho: f64,
    alpha: f64,
}

#[derive(Clone)]
enum Prop {
    Shell {
        mid: i32,
        t: f64,
        membrane: bool,
    },
    Solid {
        mid: i32,
    },
    Rod {
        mid: i32,
        area: f64,
    },
    Bar {
        mid: i32,
        area: f64,
        i1: f64,
        i2: f64,
        i12: f64,
        j: f64,
        k1: f64,
        k2: f64,
    },
}

struct BarEl {
    eid: i32,
    pid: i32,
    ga: i32,
    gb: i32,
    n1: [f64; 3],
}

enum BuiltEl {
    Std {
        eid: i32,
        kind: ElemKind,
        nodes: Vec<i32>,
        pid: i32,
    },
    Bar(BarEl),
    Conrod {
        eid: i32,
        g1: i32,
        g2: i32,
        mid: i32,
        area: f64,
    },
    Mass {
        eid: i32,
        g: i32,
        m: f64,
    },
    Spring {
        eid: i32,
        g1: i32,
        g2: i32,
        k: f64,
    },
}

enum LoadItem {
    Force {
        g: i32,
        cid: i32,
        f: [f64; 3],
        moment: bool,
    },
    Pload2 {
        p: f64,
        eids: Vec<i32>,
    },
    Pload4 {
        eid: i32,
        p: f64,
        g1: Option<i32>,
        g34: Option<i32>,
    },
    Grav {
        cid: i32,
        a: [f64; 3],
    },
    Combo {
        s: f64,
        parts: Vec<(f64, i32)>,
    },
}

struct Cord {
    o: [f64; 3],
    ex: [f64; 3],
    ey: [f64; 3],
    ez: [f64; 3],
}

struct SpcTerm {
    g: i32,
    comps: Vec<usize>,
    val: f64,
}

fn build_model(sol: i32, id_title: String, cases: &[CaseCtrl], cards: &[Vec<String>]) -> Result<Model> {
    let mut model = Model::new();
    let mut wtmass = 1.0;
    let mut nmodes = 5usize;
    let mut eig: HashMap<i32, usize> = HashMap::new();
    let mut grids: Vec<(i32, i32, [f64; 3], i32, String)> = Vec::new();
    let mut grd_cp: Option<i32> = None;
    let mut grd_cd: Option<i32> = None;
    let mut grd_ps = String::new();
    let mut cords_raw: Vec<(i32, i32, [f64; 3], [f64; 3], [f64; 3])> = Vec::new();
    let mut mats: HashMap<i32, MatRec> = HashMap::new();
    let mut props: HashMap<i32, Prop> = HashMap::new();
    let mut elements: Vec<BuiltEl> = Vec::new();
    let mut loads: HashMap<i32, Vec<LoadItem>> = HashMap::new();
    let mut spc: HashMap<i32, Vec<SpcTerm>> = HashMap::new();
    let mut spcadd: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut rbes: Vec<(i32, i32, Vec<usize>, Vec<i32>)> = Vec::new();
    let mut unknown: HashMap<String, usize> = HashMap::new();
    let mut warn_theta = false;
    let mut warn_cd = false;

    for c in cards {
        let name = c.first().map(|s| s.as_str()).unwrap_or("");
        let d = &c[1..];
        match name {
            "PARAM" => {
                let key = field(d, 0).to_ascii_uppercase();
                if key == "WTMASS" {
                    if let Some(v) = field_f64(d, 1) {
                        wtmass = v;
                    }
                }
            }
            "DEBUG" | "EIGRL" => {
                if name == "EIGRL" {
                    let sid = req_i32(d, 0, "EIGRL SID")?;
                    let nd = field_i32(d, 3).unwrap_or(5).max(1) as usize;
                    eig.insert(sid, nd);
                }
            }
            "GRDSET" => {
                grd_cp = field_i32(d, 1).or(grd_cp);
                grd_cd = field_i32(d, 5).or(grd_cd);
                if !field(d, 6).is_empty() {
                    grd_ps = field(d, 6).to_string();
                }
            }
            "GRID" => {
                let id = req_i32(d, 0, "GRID ID")?;
                let cp = field_i32(d, 1).or(grd_cp).unwrap_or(0);
                let x = [
                    field_f64(d, 2).unwrap_or(0.0),
                    field_f64(d, 3).unwrap_or(0.0),
                    field_f64(d, 4).unwrap_or(0.0),
                ];
                let cd = field_i32(d, 5).or(grd_cd).unwrap_or(0);
                let ps = if field(d, 6).is_empty() {
                    grd_ps.clone()
                } else {
                    field(d, 6).to_string()
                };
                if cd != 0 {
                    warn_cd = true;
                }
                grids.push((id, cp, x, cd, ps));
            }
            "CORD2R" => {
                let cid = req_i32(d, 0, "CORD2R")?;
                let rid = field_i32(d, 1).unwrap_or(0);
                let a = pt(d, 2)?;
                let b = pt(d, 5)?;
                let cpt = pt(d, 8)?;
                cords_raw.push((cid, rid, a, b, cpt));
            }
            "MAT1" => {
                let mid = req_i32(d, 0, "MAT1")?;
                let e = field_f64(d, 1).unwrap_or(0.0);
                let g = field_f64(d, 2);
                let nu_in = field_f64(d, 3);
                let rho = field_f64(d, 4).unwrap_or(0.0);
                let alpha = field_f64(d, 5).unwrap_or(0.0);
                let nu = match (nu_in, g) {
                    (Some(nu), _) => nu,
                    (None, Some(gv)) if e > 0.0 && gv > 0.0 => e / (2.0 * gv) - 1.0,
                    _ => 0.0,
                };
                if let (Some(nu0), Some(gv)) = (nu_in, g) {
                    if e > 0.0 {
                        let g_from_nu = e / (2.0 * (1.0 + nu0));
                        if (g_from_nu - gv).abs() > 0.01 * gv.abs().max(1.0) {
                            model.warn(format!(
                                "MAT1 {mid}: G und NU widersprechen sich — NU wird verwendet."
                            ));
                        }
                    }
                }
                mats.insert(
                    mid,
                    MatRec {
                        e,
                        nu,
                        rho,
                        alpha,
                    },
                );
            }
            "PSHELL" => {
                let pid = req_i32(d, 0, "PSHELL")?;
                let mid = field_i32(d, 1).ok_or_else(|| {
                    crate::error::FemError(format!("PSHELL {pid} ohne MID1"))
                })?;
                let t = field_f64(d, 2).unwrap_or(0.0);
                let mid2_raw = field(d, 3);
                let membrane = mid2_raw == "-1";
                let bend = field_f64(d, 4);
                if let Some(b) = bend {
                    if (b - 1.0).abs() > 1e-3 {
                        model.warn(format!(
                            "PSHELL {pid}: 12I/T³={b} wird als 1 behandelt (keine entkoppelte Biegung)."
                        ));
                    }
                }
                props.insert(pid, Prop::Shell { mid, t, membrane });
            }
            "PSOLID" => {
                let pid = req_i32(d, 0, "PSOLID")?;
                let mid = req_i32(d, 1, "PSOLID MID")?;
                props.insert(pid, Prop::Solid { mid });
            }
            "PROD" => {
                let pid = req_i32(d, 0, "PROD")?;
                let mid = req_i32(d, 1, "PROD MID")?;
                let area = field_f64(d, 2).unwrap_or(0.0);
                props.insert(pid, Prop::Rod { mid, area });
            }
            "PBAR" => {
                let pid = req_i32(d, 0, "PBAR")?;
                let mid = req_i32(d, 1, "PBAR MID")?;
                let area = field_f64(d, 2).unwrap_or(0.0);
                let i1 = field_f64(d, 3).unwrap_or(0.0);
                let i2 = field_f64(d, 4).unwrap_or(0.0);
                let j = field_f64(d, 5).unwrap_or(0.0);
                // stress-recovery continuation (8 values) then K1 K2 I12
                let k1 = field_f64(d, 16).unwrap_or(0.0);
                let k2 = field_f64(d, 17).unwrap_or(0.0);
                let i12 = field_f64(d, 18).unwrap_or(0.0);
                props.insert(
                    pid,
                    Prop::Bar {
                        mid,
                        area,
                        i1,
                        i2,
                        i12,
                        j,
                        k1,
                        k2,
                    },
                );
            }
            "PBARL" => {
                let pid = req_i32(d, 0, "PBARL")?;
                let mid = req_i32(d, 1, "PBARL MID")?;
                let typ = field(d, 3).to_ascii_uppercase();
                let dims: Vec<f64> = d.iter().skip(4).filter(|s| !s.is_empty()).filter_map(|s| parse_f64(s).ok()).collect();
                let (area, i1, i2, j) = section_library(&typ, &dims).ok_or_else(|| {
                    crate::error::FemError(format!(
                        "PBARL {pid}: Typ '{typ}' nicht unterstützt (ROD, TUBE, BAR)."
                    ))
                })?;
                props.insert(
                    pid,
                    Prop::Bar {
                        mid,
                        area,
                        i1,
                        i2,
                        i12: 0.0,
                        j,
                        k1: 0.0,
                        k2: 0.0,
                    },
                );
            }
            "CROD" => {
                elements.push(BuiltEl::Std {
                    eid: req_i32(d, 0, "CROD")?,
                    kind: ElemKind::Truss2,
                    nodes: vec![req_i32(d, 2, "CROD G1")?, req_i32(d, 3, "CROD G2")?],
                    pid: req_i32(d, 1, "CROD PID")?,
                });
            }
            "CONROD" => {
                elements.push(BuiltEl::Conrod {
                    eid: req_i32(d, 0, "CONROD")?,
                    g1: req_i32(d, 1, "CONROD G1")?,
                    g2: req_i32(d, 2, "CONROD G2")?,
                    mid: req_i32(d, 3, "CONROD MID")?,
                    area: field_f64(d, 4).unwrap_or(0.0),
                });
            }
            "CBAR" | "CBEAM" => {
                let eid = req_i32(d, 0, name)?;
                let pid = req_i32(d, 1, name)?;
                let ga = req_i32(d, 2, name)?;
                let gb = req_i32(d, 3, name)?;
                let n1 = bar_orient_fields(d);
                elements.push(BuiltEl::Bar(BarEl { eid, pid, ga, gb, n1 }));
            }
            "CQUAD4" => {
                if field_f64(d, 6).unwrap_or(0.0).abs() > 1e-8 {
                    warn_theta = true;
                }
                elements.push(BuiltEl::Std {
                    eid: req_i32(d, 0, "CQUAD4")?,
                    kind: ElemKind::Shell4,
                    nodes: vec![
                        req_i32(d, 2, "CQUAD4")?,
                        req_i32(d, 3, "CQUAD4")?,
                        req_i32(d, 4, "CQUAD4")?,
                        req_i32(d, 5, "CQUAD4")?,
                    ],
                    pid: req_i32(d, 1, "CQUAD4 PID")?,
                });
            }
            "CTRIA3" => {
                elements.push(BuiltEl::Std {
                    eid: req_i32(d, 0, "CTRIA3")?,
                    kind: ElemKind::Shell3,
                    nodes: vec![
                        req_i32(d, 2, "CTRIA3")?,
                        req_i32(d, 3, "CTRIA3")?,
                        req_i32(d, 4, "CTRIA3")?,
                    ],
                    pid: req_i32(d, 1, "CTRIA3 PID")?,
                });
            }
            "CTETRA" => {
                let eid = req_i32(d, 0, "CTETRA")?;
                let pid = req_i32(d, 1, "CTETRA PID")?;
                let nodes = ints_skip(d, 2, 10);
                let kind = if nodes.len() >= 10 {
                    ElemKind::Tet10
                } else {
                    ElemKind::Tet4
                };
                elements.push(BuiltEl::Std {
                    eid,
                    kind,
                    nodes: nodes[..kind.nnodes()].to_vec(),
                    pid,
                });
            }
            "CHEXA" => {
                let eid = req_i32(d, 0, "CHEXA")?;
                let pid = req_i32(d, 1, "CHEXA PID")?;
                let nodes = ints_skip(d, 2, 20);
                let kind = if nodes.len() >= 20 {
                    ElemKind::Hex20
                } else {
                    ElemKind::Hex8
                };
                elements.push(BuiltEl::Std {
                    eid,
                    kind,
                    nodes: nodes[..kind.nnodes()].to_vec(),
                    pid,
                });
            }
            "CPENTA" => {
                let eid = req_i32(d, 0, "CPENTA")?;
                let pid = req_i32(d, 1, "CPENTA PID")?;
                let nodes = ints_skip(d, 2, 15);
                let kind = if nodes.len() >= 15 {
                    ElemKind::Wedge15
                } else {
                    ElemKind::Wedge6
                };
                elements.push(BuiltEl::Std {
                    eid,
                    kind,
                    nodes: nodes[..kind.nnodes()].to_vec(),
                    pid,
                });
            }
            "CELAS2" => {
                let eid = req_i32(d, 0, "CELAS2")?;
                let k = field_f64(d, 1).unwrap_or(0.0);
                let g1 = field_i32(d, 2).unwrap_or(0);
                let g2 = field_i32(d, 4).unwrap_or(0);
                if g1 == 0 || g2 == 0 {
                    model.warn(format!("CELAS2 {eid}: Feder gegen Erde wird nicht unterstützt."));
                } else {
                    elements.push(BuiltEl::Spring { eid, g1, g2, k });
                }
            }
            "CONM2" => {
                let eid = req_i32(d, 0, "CONM2")?;
                let g = req_i32(d, 1, "CONM2 G")?;
                let m = field_f64(d, 3).unwrap_or(0.0);
                let off = field_f64(d, 4).unwrap_or(0.0).hypot(field_f64(d, 5).unwrap_or(0.0)).hypot(field_f64(d, 6).unwrap_or(0.0));
                if off > 1e-12 {
                    model.warn(format!("CONM2 {eid}: Offset wird ignoriert."));
                }
                elements.push(BuiltEl::Mass { eid, g, m });
            }
            "RBE2" => {
                let eid = req_i32(d, 0, "RBE2")?;
                let gn = req_i32(d, 1, "RBE2 GN")?;
                let cm = comps(field(d, 2));
                let dep = expand_ids(&d[3..]);
                rbes.push((eid, gn, cm, dep));
            }
            "FORCE" | "MOMENT" => {
                let sid = req_i32(d, 0, name)?;
                let g = req_i32(d, 1, name)?;
                let cid = field_i32(d, 2).unwrap_or(0);
                let mag = field_f64(d, 3).unwrap_or(0.0);
                let n = [
                    field_f64(d, 4).unwrap_or(0.0),
                    field_f64(d, 5).unwrap_or(0.0),
                    field_f64(d, 6).unwrap_or(0.0),
                ];
                let ln = norm3(n);
                let f = if ln < 1e-30 {
                    [0.0; 3]
                } else {
                    scale3(n, mag / ln)
                };
                loads.entry(sid).or_default().push(LoadItem::Force {
                    g,
                    cid,
                    f,
                    moment: name == "MOMENT",
                });
            }
            "PLOAD2" => {
                let sid = req_i32(d, 0, "PLOAD2")?;
                let p = field_f64(d, 1).unwrap_or(0.0);
                let eids = expand_ids(&d[2..]);
                loads.entry(sid).or_default().push(LoadItem::Pload2 { p, eids });
            }
            "PLOAD4" => {
                let sid = req_i32(d, 0, "PLOAD4")?;
                let eid = req_i32(d, 1, "PLOAD4 EID")?;
                let p1 = field_f64(d, 2).unwrap_or(0.0);
                let p = match (field_f64(d, 3), field_f64(d, 4), field_f64(d, 5)) {
                    (None, None, None) => p1,
                    (a, b, c) => 0.25 * (p1 + a.unwrap_or(p1) + b.unwrap_or(p1) + c.unwrap_or(p1)),
                };
                let g1 = field_i32(d, 6);
                let g34 = field_i32(d, 7);
                loads.entry(sid).or_default().push(LoadItem::Pload4 { eid, p, g1, g34 });
            }
            "GRAV" => {
                let sid = req_i32(d, 0, "GRAV")?;
                let cid = field_i32(d, 1).unwrap_or(0);
                let a = field_f64(d, 2).unwrap_or(0.0);
                let n = [
                    field_f64(d, 3).unwrap_or(0.0),
                    field_f64(d, 4).unwrap_or(0.0),
                    field_f64(d, 5).unwrap_or(0.0),
                ];
                let ln = norm3(n);
                let dir = if ln < 1e-30 { [0.0, 0.0, -1.0] } else { scale3(n, 1.0 / ln) };
                loads.entry(sid).or_default().push(LoadItem::Grav {
                    cid,
                    a: scale3(dir, a),
                });
            }
            "LOAD" => {
                let sid = req_i32(d, 0, "LOAD")?;
                let s = field_f64(d, 1).unwrap_or(1.0);
                let mut parts = Vec::new();
                let mut i = 2;
                while i + 1 < d.len() {
                    if d[i].is_empty() && d[i + 1].is_empty() {
                        i += 2;
                        continue;
                    }
                    let si = parse_f64(&d[i]).unwrap_or(0.0);
                    let li = parse_i32(&d[i + 1]).unwrap_or(0);
                    if li != 0 {
                        parts.push((si, li));
                    }
                    i += 2;
                }
                loads.entry(sid).or_default().push(LoadItem::Combo { s, parts });
            }
            "SPC" => {
                let sid = req_i32(d, 0, "SPC")?;
                let g = req_i32(d, 1, "SPC G")?;
                let c = comps(field(d, 2));
                let val = field_f64(d, 3).unwrap_or(0.0);
                spc.entry(sid).or_default().push(SpcTerm { g, comps: c, val });
            }
            "SPC1" => {
                let sid = req_i32(d, 0, "SPC1")?;
                let c = comps(field(d, 1));
                let ids = expand_ids(&d[2..]);
                let entry = spc.entry(sid).or_default();
                for g in ids {
                    entry.push(SpcTerm {
                        g,
                        comps: c.clone(),
                        val: 0.0,
                    });
                }
            }
            "SPCADD" => {
                let sid = req_i32(d, 0, "SPCADD")?;
                let ids = d.iter().filter(|s| !s.is_empty()).filter_map(|s| parse_i32(s).ok()).collect();
                spcadd.insert(sid, ids);
            }
            "" => {}
            other => {
                if let Some(card) = mystran_manifest::lookup(other) {
                    if card.status == CardStatus::Declined {
                        return err(format!("{other}: {}", card.note));
                    }
                }
                *unknown.entry(other.to_string()).or_insert(0) += 1;
            }
        }
    }

    if warn_theta {
        model.warn("CQUAD4-Materialwinkel wird ignoriert (isotropes MITC4).");
    }
    if warn_cd {
        model.warn("GRID CD≠0 wird ignoriert — Freiheitsgrade bleiben im Basissystem.");
    }
    if !unknown.is_empty() {
        let mut names: Vec<_> = unknown.keys().cloned().collect();
        names.sort();
        model.warn(format!(
            "Nicht unterstützte Bulk-Karten ignoriert: {}",
            names.join(", ")
        ));
    }

    let cords = resolve_cords(&cords_raw)?;
    let mut seen_n = HashSet::new();
    for (id, cp, x, _cd, ps) in &grids {
        if !seen_n.insert(*id) {
            return err(format!("GRID {id} doppelt."));
        }
        let xb = if *cp == 0 {
            *x
        } else {
            let c = cords.get(cp).ok_or_else(|| {
                crate::error::FemError(format!("GRID {id}: Koordinatensystem {cp} fehlt."))
            })?;
            [
                c.o[0] + x[0] * c.ex[0] + x[1] * c.ey[0] + x[2] * c.ez[0],
                c.o[1] + x[0] * c.ex[1] + x[1] * c.ey[1] + x[2] * c.ez[1],
                c.o[2] + x[0] * c.ex[2] + x[1] * c.ey[2] + x[2] * c.ez[2],
            ]
        };
        model.id_to_index.insert(*id, model.coords.len());
        model.node_ids.push(*id);
        model.coords.push(xb);
        let _ = ps;
    }

    let mut elem_kind: HashMap<i32, ElemKind> = HashMap::new();
    let mut elem_nodes: HashMap<i32, Vec<i32>> = HashMap::new();

    for el in &elements {
        match el {
            BuiltEl::Std { eid, kind, nodes, pid } => {
                let (elset, membrane) = bind_prop(&mut model, &mats, &props, *pid, wtmass)?;
                let kind = if membrane {
                    match kind {
                        ElemKind::Shell4 => ElemKind::Mem4,
                        ElemKind::Shell3 => ElemKind::Mem3,
                        k => *k,
                    }
                } else {
                    *kind
                };
                model.elements.push(Element {
                    id: *eid,
                    kind,
                    nodes: nodes.clone(),
                    elset,
                });
                elem_kind.insert(*eid, kind);
                elem_nodes.insert(*eid, nodes.clone());
            }
            BuiltEl::Conrod { eid, g1, g2, mid, area } => {
                let elset = format!("CONROD{eid}");
                bind_mat(&mut model, &mats, &elset, *mid, wtmass)?;
                model.elset_thickness.insert(elset.clone(), *area);
                model.elements.push(Element {
                    id: *eid,
                    kind: ElemKind::Truss2,
                    nodes: vec![*g1, *g2],
                    elset,
                });
                elem_kind.insert(*eid, ElemKind::Truss2);
                elem_nodes.insert(*eid, vec![*g1, *g2]);
            }
            BuiltEl::Bar(b) => {
                let prop = props.get(&b.pid).ok_or_else(|| {
                    crate::error::FemError(format!("CBAR {}: PID {} fehlt.", b.eid, b.pid))
                })?;
                let Prop::Bar { mid, area, i1, i2, i12, j, k1, k2 } = prop.clone() else {
                    return err(format!("CBAR {}: PID {} ist kein PBAR/PBARL.", b.eid, b.pid));
                };
                let elset = format!("PID{}E{}", b.pid, b.eid);
                bind_mat(&mut model, &mats, &elset, mid, wtmass)?;
                let xa = grid_xyz(&model.node_ids, &model.coords, b.ga)?;
                let xb = grid_xyz(&model.node_ids, &model.coords, b.gb)?;
                let n1 = if g0_flag(&b.n1) {
                    [0.0, 1.0, 0.0]
                } else {
                    beam_n1(xa, xb, b.n1)
                };
                let mut sec = BeamSection::general(area, i2, i12, i1, j, n1);
                sec.k11 = if k1 > 0.0 { k1 } else { 1.0e6 };
                sec.k22 = if k2 > 0.0 { k2 } else { 1.0e6 };
                model.elset_beam.insert(elset.clone(), sec);
                model.elements.push(Element {
                    id: b.eid,
                    kind: ElemKind::Beam31,
                    nodes: vec![b.ga, b.gb],
                    elset,
                });
                elem_kind.insert(b.eid, ElemKind::Beam31);
                elem_nodes.insert(b.eid, vec![b.ga, b.gb]);
            }
            BuiltEl::Mass { eid, g, m } => {
                let elset = format!("MASS{eid}");
                model.elset_mass.insert(elset.clone(), *m);
                model.elements.push(Element {
                    id: *eid,
                    kind: ElemKind::Mass,
                    nodes: vec![*g],
                    elset,
                });
            }
            BuiltEl::Spring { eid, g1, g2, k } => {
                let elset = format!("SPR{eid}");
                model.elset_spring.insert(elset.clone(), *k);
                model.elements.push(Element {
                    id: *eid,
                    kind: ElemKind::SpringA,
                    nodes: vec![*g1, *g2],
                    elset,
                });
            }
        }
    }

    // G0 orientation: re-read bars whose n1.x is a grid id encoded as NaN-free sentinel.
    // Implemented in orient_g0 pass below if we stored g0 in n1[0] and n1[1] as a flag.
    apply_g0(&mut model, &elements)?;

    let mut slave_nodes = HashSet::new();
    let mut ref_nodes = HashSet::new();
    for (eid, gn, cm, dep) in &rbes {
        if !cm.contains(&1) || !cm.contains(&2) || !cm.contains(&3) {
            model.warn(format!("RBE2 {eid}: nur CM mit 123 (und optional 456) wird als Starrkörper übernommen."));
        }
        let name = format!("RBE{eid}");
        model.nsets.insert(name.clone(), dep.clone());
        model.rigid_bodies.push(RigidBody {
            nset: name,
            ref_node: *gn,
            rot_node: None,
        });
        ref_nodes.insert(*gn);
        slave_nodes.extend(dep.iter().copied());
    }

    let pin_rot: Vec<i32> = if rbes.is_empty() {
        Vec::new()
    } else {
        let structural: HashSet<i32> = model
            .elements
            .iter()
            .filter(|e| e.kind.is_shell() || e.kind.is_beam())
            .flat_map(|e| e.nodes.iter().copied())
            .collect();
        model
            .node_ids
            .iter()
            .copied()
            .filter(|id| !ref_nodes.contains(id) && !structural.contains(id) && !slave_nodes.contains(id))
            .collect()
    };

    let heading = cases
        .first()
        .map(|c| {
            if !c.title.is_empty() {
                c.title.clone()
            } else if !id_title.is_empty() {
                id_title.clone()
            } else {
                "MYSTRAN".into()
            }
        })
        .unwrap_or_else(|| "MYSTRAN".into());
    model.heading = heading;

    let procedure = match sol {
        3 | 103 => {
            let nd = cases
                .first()
                .and_then(|c| c.method)
                .and_then(|m| eig.get(&m).copied())
                .unwrap_or(nmodes);
            nmodes = nd;
            Procedure::Frequency { nmodes }
        }
        1 | 101 | _ => {
            if !matches!(sol, 1 | 101) {
                model.warn(format!("SOL {sol} wird als lineare Statik (SOL 101) gerechnet."));
            }
            Procedure::Static {
                nlgeom: false,
                increments: 1,
                riks: false,
            }
        }
    };

    let permanent: Vec<(i32, String)> = grids.iter().map(|(id, _, _, _, ps)| (*id, ps.clone())).collect();

    for (ci, case) in cases.iter().enumerate() {
        let bc_from = model.bcs.len();
        let c_from = model.cloads.len();
        let d_from = model.dloads.len();
        let mut have: HashSet<(i32, usize)> = HashSet::new();
        for (id, ps) in &permanent {
            push_spc(&mut model, &mut have, *id, &comps(ps), 0.0);
        }
        for id in &pin_rot {
            push_spc(&mut model, &mut have, *id, &[4, 5, 6], 0.0);
        }
        if let Some(sid) = case.spc {
            let terms = resolve_spc(sid, &spc, &spcadd)?;
            for t in terms {
                push_spc(&mut model, &mut have, t.g, &t.comps, t.val);
            }
        }
        if let Some(sid) = case.load {
            apply_load(
                &mut model,
                sid,
                1.0,
                &loads,
                &cords,
                &elem_kind,
                &elem_nodes,
                &mut HashSet::new(),
            )?;
        }
        let label = if !case.label.is_empty() {
            case.label.clone()
        } else if !case.subtitle.is_empty() {
            case.subtitle.clone()
        } else if cases.len() > 1 {
            format!("SUBCASE {}", ci + 1)
        } else {
            case.title.clone()
        };
        model.case_labels.push(label);
        model.steps.push(AnalysisStep {
            procedure: procedure.clone(),
            n_cload: model.cloads.len(),
            n_dload: model.dloads.len(),
            n_bc: model.bcs.len(),
            cload_from: c_from,
            dload_from: d_from,
            bc_from,
        });
    }
    if model.steps.len() > 1 {
        model.independent_steps = true;
    }
    model.procedure = procedure;

    if model.node_ids.is_empty() {
        return err("Keine GRID-Karten im MYSTRAN-Deck.");
    }
    if model.elements.is_empty() {
        return err("Keine Elemente im MYSTRAN-Deck.");
    }
    model.compact();
    let _ = nmodes;
    Ok(model)
}

fn bind_prop(
    model: &mut Model,
    mats: &HashMap<i32, MatRec>,
    props: &HashMap<i32, Prop>,
    pid: i32,
    wtmass: f64,
) -> Result<(String, bool)> {
    let prop = props.get(&pid).ok_or_else(|| {
        crate::error::FemError(format!("Eigenschaft {pid} fehlt."))
    })?;
    let elset = format!("PID{pid}");
    let membrane = match prop {
        Prop::Shell { mid, t, membrane } => {
            bind_mat(model, mats, &elset, *mid, wtmass)?;
            model.elset_thickness.insert(elset.clone(), *t);
            *membrane
        }
        Prop::Solid { mid } => {
            bind_mat(model, mats, &elset, *mid, wtmass)?;
            false
        }
        Prop::Rod { mid, area } => {
            bind_mat(model, mats, &elset, *mid, wtmass)?;
            model.elset_thickness.insert(elset.clone(), *area);
            false
        }
        Prop::Bar { .. } => false,
    };
    Ok((elset, membrane))
}

fn bind_mat(model: &mut Model, mats: &HashMap<i32, MatRec>, elset: &str, mid: i32, wtmass: f64) -> Result<()> {
    let m = mats.get(&mid).ok_or_else(|| {
        crate::error::FemError(format!("MAT1 {mid} fehlt."))
    })?;
    let name = format!("M{mid}");
    model.materials.entry(name.clone()).or_insert(Material {
        e: m.e,
        nu: m.nu,
        density: m.rho * wtmass,
        alpha: m.alpha,
        ..Material::default()
    });
    model.elset_material.insert(elset.to_string(), name);
    Ok(())
}

fn grid_xyz(ids: &[i32], coords: &[[f64; 3]], id: i32) -> Result<[f64; 3]> {
    let i = ids
        .iter()
        .position(|&n| n == id)
        .ok_or_else(|| crate::error::FemError(format!("Knoten {id} fehlt.")))?;
    Ok(coords[i])
}

fn apply_g0(model: &mut Model, elements: &[BuiltEl]) -> Result<()> {
    for el in elements {
        let BuiltEl::Bar(b) = el else { continue };
        if !g0_flag(&b.n1) {
            continue;
        }
        let g0 = b.n1[0] as i32;
        let ia = model.node_index(b.ga)?;
        let ib = model.node_index(b.gb)?;
        let ig = model.node_index(g0)?;
        let xa = model.coords[ia];
        let xb = model.coords[ib];
        let xg = model.coords[ig];
        let v = [xg[0] - xa[0], xg[1] - xa[1], xg[2] - xa[2]];
        let n1 = beam_n1(xa, xb, v);
        let elset = format!("PID{}E{}", b.pid, b.eid);
        if let Some(sec) = model.elset_beam.get_mut(&elset) {
            sec.n1 = n1;
        }
    }
    Ok(())
}

fn g0_flag(n1: &[f64; 3]) -> bool {
    n1[2].is_nan()
}

fn bar_orient_fields(d: &[String]) -> [f64; 3] {
    let f5 = field(d, 4);
    let f6 = field(d, 5);
    let f7 = field(d, 6);
    if !f5.is_empty() && f6.is_empty() && f7.is_empty() {
        if let Ok(g0) = parse_i32(f5) {
            return [g0 as f64, 0.0, f64::NAN];
        }
    }
    [
        field_f64(d, 4).unwrap_or(0.0),
        field_f64(d, 5).unwrap_or(0.0),
        field_f64(d, 6).unwrap_or(0.0),
    ]
}

fn beam_n1(xa: [f64; 3], xb: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    let t = sub3(xb, xa);
    let tn = norm3(t).max(1e-30);
    let t = scale3(t, 1.0 / tn);
    let mut y = sub3(v, scale3(t, dot3(v, t)));
    if norm3(y) < 1e-12 {
        let alt = if t[2].abs() < 0.9 { [0.0, 0.0, 1.0] } else { [0.0, 1.0, 0.0] };
        y = sub3(alt, scale3(t, dot3(alt, t)));
    }
    let yn = norm3(y).max(1e-30);
    scale3(y, 1.0 / yn)
}

fn section_library(typ: &str, dims: &[f64]) -> Option<(f64, f64, f64, f64)> {
    let pi = std::f64::consts::PI;
    match typ {
        "ROD" => {
            let r = *dims.first()?;
            let a = pi * r * r;
            let i = pi * r.powi(4) / 4.0;
            Some((a, i, i, 2.0 * i))
        }
        "TUBE" => {
            let d1 = *dims.first()?;
            let d2 = *dims.get(1)?;
            let (ro, ri) = if d1 >= d2 { (d1, d2) } else { (d2, d1) };
            let a = pi * (ro * ro - ri * ri);
            let i = pi / 4.0 * (ro.powi(4) - ri.powi(4));
            Some((a, i, i, 2.0 * i))
        }
        "BAR" => {
            let w = *dims.first()?;
            let h = *dims.get(1)?;
            let sec = BeamSection::rect(w, h, [0.0, 1.0, 0.0]);
            // i22 of rect is Nastran I1 (bending in the y-plane), i11 is I2.
            Some((sec.area, sec.i22, sec.i11, sec.jtor))
        }
        _ => None,
    }
}

fn resolve_cords(raw: &[(i32, i32, [f64; 3], [f64; 3], [f64; 3])]) -> Result<HashMap<i32, Cord>> {
    let mut out = HashMap::new();
    let mut left: Vec<_> = raw.to_vec();
    while !left.is_empty() {
        let n0 = left.len();
        left.retain(|(cid, rid, a, b, c)| {
            if *rid != 0 && !out.contains_key(rid) {
                return true;
            }
            let (a, b, c) = if *rid == 0 {
                (*a, *b, *c)
            } else {
                let p = &out[rid];
                (to_basic(p, *a), to_basic(p, *b), to_basic(p, *c))
            };
            let ez = unit(sub3(b, a));
            let ac = sub3(c, a);
            let ex = unit(sub3(ac, scale3(ez, dot3(ac, ez))));
            let ey = cross3(ez, ex);
            out.insert(
                *cid,
                Cord {
                    o: a,
                    ex,
                    ey,
                    ez,
                },
            );
            false
        });
        if left.len() == n0 {
            return err(format!(
                "CORD2R {} hängt an einem fehlenden System.",
                left[0].0
            ));
        }
    }
    Ok(out)
}

fn to_basic(c: &Cord, x: [f64; 3]) -> [f64; 3] {
    [
        c.o[0] + x[0] * c.ex[0] + x[1] * c.ey[0] + x[2] * c.ez[0],
        c.o[1] + x[0] * c.ex[1] + x[1] * c.ey[1] + x[2] * c.ez[1],
        c.o[2] + x[0] * c.ex[2] + x[1] * c.ey[2] + x[2] * c.ez[2],
    ]
}

fn resolve_spc<'a>(
    sid: i32,
    spc: &'a HashMap<i32, Vec<SpcTerm>>,
    add: &HashMap<i32, Vec<i32>>,
) -> Result<Vec<&'a SpcTerm>> {
    fn walk<'a>(
        sid: i32,
        spc: &'a HashMap<i32, Vec<SpcTerm>>,
        add: &HashMap<i32, Vec<i32>>,
        seen: &mut HashSet<i32>,
        out: &mut Vec<&'a SpcTerm>,
    ) -> Result<()> {
        if !seen.insert(sid) {
            return err(format!("SPCADD {sid} ist zyklisch."));
        }
        if let Some(ids) = add.get(&sid) {
            for s in ids {
                walk(*s, spc, add, seen, out)?;
            }
        }
        if let Some(terms) = spc.get(&sid) {
            out.extend(terms.iter());
        }
        if !add.contains_key(&sid) && !spc.contains_key(&sid) {
            return err(format!("SPC-Set {sid} fehlt."));
        }
        Ok(())
    }
    let mut out = Vec::new();
    walk(sid, spc, add, &mut HashSet::new(), &mut out)?;
    Ok(out)
}

fn apply_load(
    model: &mut Model,
    sid: i32,
    scale: f64,
    loads: &HashMap<i32, Vec<LoadItem>>,
    cords: &HashMap<i32, Cord>,
    elem_kind: &HashMap<i32, ElemKind>,
    elem_nodes: &HashMap<i32, Vec<i32>>,
    seen: &mut HashSet<i32>,
) -> Result<()> {
    if !seen.insert(sid) {
        return err(format!("LOAD {sid} ist zyklisch."));
    }
    let items = loads.get(&sid).ok_or_else(|| {
        crate::error::FemError(format!("Lastsatz {sid} fehlt."))
    })?;
    for it in items {
        match it {
            LoadItem::Combo { s, parts } => {
                for (si, li) in parts {
                    apply_load(model, *li, scale * s * si, loads, cords, elem_kind, elem_nodes, seen)?;
                }
            }
            LoadItem::Force { g, cid, f, moment } => {
                let v = scale3(vec_in_basic(*cid, *f, cords)?, scale);
                let base = if *moment { 3 } else { 0 };
                for k in 0..3 {
                    if v[k].abs() > 0.0 {
                        model.cloads.push(Cload {
                            node: *g,
                            dof: base + k,
                            mag: v[k],
                            amplitude: String::new(),
                        });
                    }
                }
            }
            LoadItem::Grav { cid, a } => {
                let v = scale3(vec_in_basic(*cid, *a, cords)?, scale);
                let mag = norm3(v);
                if mag > 0.0 {
                    model.dloads.push(Dload::Grav { mag, dir: v });
                }
            }
            LoadItem::Pload2 { p, eids } => {
                for eid in eids {
                    push_pressure(model, *eid, -scale * *p, None, None, elem_kind, elem_nodes)?;
                }
            }
            LoadItem::Pload4 { eid, p, g1, g34 } => {
                let kind = elem_kind.get(eid).copied();
                let shell = kind.map(|k| k.is_shell() || k.is_membrane()).unwrap_or(false);
                let mag = if shell { -scale * *p } else { scale * *p };
                push_pressure(model, *eid, mag, *g1, *g34, elem_kind, elem_nodes)?;
            }
        }
    }
    seen.remove(&sid);
    Ok(())
}

fn push_pressure(
    model: &mut Model,
    eid: i32,
    mag: f64,
    g1: Option<i32>,
    g34: Option<i32>,
    elem_kind: &HashMap<i32, ElemKind>,
    elem_nodes: &HashMap<i32, Vec<i32>>,
) -> Result<()> {
    let kind = elem_kind.get(&eid).copied().ok_or_else(|| {
        crate::error::FemError(format!("Druck auf unbekanntes Element {eid}."))
    })?;
    let face = if kind.is_shell() || kind.is_membrane() || kind.is_truss() {
        1
    } else if let (Some(a), Some(b)) = (g1, g34) {
        let nodes = elem_nodes.get(&eid).ok_or_else(|| {
            crate::error::FemError(format!("Element {eid} ohne Knoten."))
        })?;
        solid_face(kind, nodes, a, b).ok_or_else(|| {
            crate::error::FemError(format!("PLOAD4 {eid}: Fläche {a}/{b} nicht gefunden."))
        })?
    } else if kind == ElemKind::Hex8 || kind == ElemKind::Tet4 || kind == ElemKind::Wedge6 {
        return err(format!(
            "PLOAD4 auf Volumenelement {eid} braucht G1 und G34 zur Flächenwahl."
        ));
    } else {
        1
    };
    model.dloads.push(Dload::Pressure { elem: eid, face, mag });
    Ok(())
}

fn solid_face(kind: ElemKind, nodes: &[i32], g1: i32, g34: i32) -> Option<i32> {
    let faces: &[&[usize]] = match kind {
        ElemKind::Hex8 | ElemKind::Hex20 => &[
            &[0, 1, 2, 3],
            &[4, 5, 6, 7],
            &[0, 1, 5, 4],
            &[1, 2, 6, 5],
            &[2, 3, 7, 6],
            &[3, 0, 4, 7],
        ],
        ElemKind::Tet4 | ElemKind::Tet10 => &[&[0, 1, 2], &[0, 3, 1], &[1, 3, 2], &[2, 3, 0]],
        ElemKind::Wedge6 | ElemKind::Wedge15 => &[
            &[0, 1, 2],
            &[3, 4, 5],
            &[0, 1, 4, 3],
            &[1, 2, 5, 4],
            &[2, 0, 3, 5],
        ],
        _ => return None,
    };
    // CalculiX hex faces are not the same order as this list for P3.. — match by membership only.
    let hex_ccx: Option<[[usize; 4]; 6]> = if matches!(kind, ElemKind::Hex8 | ElemKind::Hex20) {
        Some([
            [0, 1, 2, 3],
            [4, 7, 6, 5],
            [0, 4, 5, 1],
            [1, 5, 6, 2],
            [2, 6, 7, 3],
            [3, 7, 4, 0],
        ])
    } else {
        None
    };
    if let Some(ff) = hex_ccx {
        for (i, f) in ff.iter().enumerate() {
            let ids: Vec<i32> = f.iter().filter_map(|k| nodes.get(*k).copied()).collect();
            if ids.contains(&g1) && ids.contains(&g34) {
                return Some(i as i32 + 1);
            }
        }
        return None;
    }
    let _ = faces;
    match kind {
        ElemKind::Tet4 | ElemKind::Tet10 => {
            let ff: [[usize; 3]; 4] = [[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]];
            for (i, f) in ff.iter().enumerate() {
                let ids: Vec<i32> = f.iter().filter_map(|k| nodes.get(*k).copied()).collect();
                if ids.contains(&g1) && ids.contains(&g34) {
                    return Some((i as i32) + 1);
                }
            }
            None
        }
        _ => None,
    }
}

fn vec_in_basic(cid: i32, v: [f64; 3], cords: &HashMap<i32, Cord>) -> Result<[f64; 3]> {
    if cid == 0 {
        return Ok(v);
    }
    let c = cords.get(&cid).ok_or_else(|| {
        crate::error::FemError(format!("Koordinatensystem {cid} fehlt."))
    })?;
    Ok([
        v[0] * c.ex[0] + v[1] * c.ey[0] + v[2] * c.ez[0],
        v[0] * c.ex[1] + v[1] * c.ey[1] + v[2] * c.ez[1],
        v[0] * c.ex[2] + v[1] * c.ey[2] + v[2] * c.ez[2],
    ])
}

fn push_spc(model: &mut Model, have: &mut HashSet<(i32, usize)>, node: i32, comps: &[usize], val: f64) {
    for &c in comps {
        if !(1..=6).contains(&c) {
            continue;
        }
        let dof = c - 1;
        if have.insert((node, dof)) {
            model.bcs.push(Boundary { node, dof, value: val });
        }
    }
}

fn comps(s: &str) -> Vec<usize> {
    let t = s.trim();
    if t.is_empty() {
        return Vec::new();
    }
    let digits = if let Ok(v) = parse_i32(t) {
        v.abs().to_string()
    } else {
        t.to_string()
    };
    digits
        .chars()
        .filter_map(|c| c.to_digit(10))
        .map(|d| d as usize)
        .filter(|d| (1..=6).contains(d))
        .collect()
}

fn expand_ids(fields: &[String]) -> Vec<i32> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < fields.len() {
        if fields[i].is_empty() {
            i += 1;
            continue;
        }
        if fields[i].eq_ignore_ascii_case("THRU") {
            i += 1;
            continue;
        }
        let Ok(id) = parse_i32(&fields[i]) else {
            i += 1;
            continue;
        };
        if i + 2 <= fields.len() - 1 && fields.get(i + 1).map(|s| s.eq_ignore_ascii_case("THRU")).unwrap_or(false)
        {
            if let Ok(b) = parse_i32(&fields[i + 2]) {
                let (lo, hi) = if id <= b { (id, b) } else { (b, id) };
                for g in lo..=hi {
                    out.push(g);
                }
                i += 3;
                continue;
            }
        }
        out.push(id);
        i += 1;
    }
    out
}

fn ints_skip(d: &[String], from: usize, n: usize) -> Vec<i32> {
    d.iter()
        .skip(from)
        .filter(|s| !s.is_empty())
        .filter_map(|s| parse_i32(s).ok())
        .take(n)
        .collect()
}

fn field(d: &[String], i: usize) -> &str {
    d.get(i).map(|s| s.as_str()).unwrap_or("")
}

fn field_f64(d: &[String], i: usize) -> Option<f64> {
    let s = field(d, i);
    if s.is_empty() {
        None
    } else {
        parse_f64(s).ok()
    }
}

fn field_i32(d: &[String], i: usize) -> Option<i32> {
    let s = field(d, i);
    if s.is_empty() {
        None
    } else {
        parse_i32(s).ok()
    }
}

fn req_i32(d: &[String], i: usize, what: &str) -> Result<i32> {
    field_i32(d, i).ok_or_else(|| crate::error::FemError(format!("{what}: Feld {} fehlt.", i + 1)))
}

fn pt(d: &[String], i: usize) -> Result<[f64; 3]> {
    Ok([
        field_f64(d, i).unwrap_or(0.0),
        field_f64(d, i + 1).unwrap_or(0.0),
        field_f64(d, i + 2).unwrap_or(0.0),
    ])
}

fn strip_dollar(line: &str) -> &str {
    match line.find('$') {
        Some(i) => &line[..i],
        None => line,
    }
}

pub fn parse_f64(s: &str) -> Result<f64> {
    let t = s.trim();
    if t.is_empty() {
        return err("leere Zahl");
    }
    let mut u = t.replace(['d', 'D'], "E");
    if !u.contains('E') && !u.contains('e') {
        // Nastran 1.0+3 / 1.0-3, but not a leading sign.
        if let Some(k) = u.rfind(['+', '-']) {
            if k > 0 && u.as_bytes()[k - 1].is_ascii_digit() {
                u.insert(k, 'E');
            }
        }
    }
    u.parse::<f64>()
        .map_err(|_| crate::error::FemError(format!("keine Zahl: {s}")))
}

fn parse_i32(s: &str) -> Result<i32> {
    let t = s.trim();
    if let Ok(v) = t.parse::<i32>() {
        return Ok(v);
    }
    let f = parse_f64(t)?;
    if (f - f.round()).abs() > 1e-6 {
        return err(format!("keine Ganzzahl: {s}"));
    }
    Ok(f.round() as i32)
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
fn unit(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a).max(1e-30);
    scale3(a, 1.0 / n)
}
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::solve;

    #[test]
    fn nastran_exponent_and_detect() {
        assert!((parse_f64("1.0+3").unwrap() - 1000.0).abs() < 1e-9);
        assert!((parse_f64("2.5-4").unwrap() - 2.5e-4).abs() < 1e-12);
        assert!((parse_f64("1.E+7").unwrap() - 1e7).abs() < 1.0);
        let deck = "SOL 1\nCEND\nBEGIN BULK\nGRID,1,,0,0,0\nENDDATA\n";
        assert!(is_mystran_deck(deck));
        assert!(!is_mystran_deck("*NODE\n1, 0, 0, 0\n"));
    }

    #[test]
    fn crod_tip_displacement() {
        let deck = r#"
ID rod
SOL 101
CEND
TITLE = rod
SPC = 1
LOAD = 1
BEGIN BULK
GRID,1,,0.,0.,0.
GRID,2,,10.,0.,0.
CROD,1,1,1,2
PROD,1,1,2.0
MAT1,1,210000.,,0.3
SPC1,1,123456,1
FORCE,1,2,,1000.,1.,0.,0.
ENDDATA
"#;
        let m = parse_with_base(deck, None).unwrap();
        assert_eq!(m.elements.len(), 1);
        assert_eq!(m.elements[0].kind, ElemKind::Truss2);
        let out = solve(m).unwrap();
        let ux = out.u[out.model.node_index(2).unwrap()][0];
        let expect = 1000.0 * 10.0 / (210000.0 * 2.0);
        assert!((ux - expect).abs() < 1e-8, "ux={ux} expect={expect}");
        let f06 = crate::f06::write_f06(&out.model, &out);
        assert!(f06.contains("D I S P L A C E M E N T S"));
        assert!(f06.contains("S P C   F O R C E S"));
        assert!(f06.contains("2.380952E-02"), "{f06}");
    }

    #[test]
    fn fixed_field_columns_and_bar_load_scale() {
        let grid = format!(
            "{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}",
            "GRID", "11", "", "0.0", "0.0", "0.0", "", "123456"
        );
        let (f, large) = split_card(&grid).unwrap();
        assert!(!large);
        assert_eq!(f[1], "11");
        assert_eq!(f[7], "123456");

        let deck = r#"
SOL 1
CEND
SPC = 1
LOAD = 1
BEGIN BULK
GRID,11,,0.,0.,0.
GRID,12,,100.,0.,0.
CBAR,1,10,11,12,0.,1.,0.
PBAR,10,1,1.0,1.0+9,20.,1.0
MAT1,1,1000.,,0.0
SPC,1,12,123456,0.
FORCE,2,11,,1.,0.,0.,1.
LOAD,1,1.0,5.0,2
ENDDATA
"#;
        let m = parse_with_base(deck, None).unwrap();
        let out = solve(m).unwrap();
        let i = out.model.node_index(11).unwrap();
        // B31 is a one-point Timoshenko beam. With a huge shear factor the
        // tip stiffness is 4EI/L³ (not the cubic 3EI/L³). I2=20 bends about
        // local y, so Fz moves uz; I1=1e9 keeps uy at rest.
        let uz = out.u[i][2];
        let uy = out.u[i][1];
        let expect = 5.0 * 100f64.powi(3) / (4.0 * 1000.0 * 20.0);
        assert!(uy.abs() < 1e-6, "uy={uy}");
        assert!(
            (uz - expect).abs() / expect < 1e-6,
            "uz={uz} expect={expect}"
        );
    }

    #[test]
    fn pload2_opposes_normal_and_two_subcases() {
        let deck = r#"
SOL 1
CEND
SPC = 1
SUBCASE 1
  LOAD = 1
SUBCASE 2
  LOAD = 2
BEGIN BULK
GRID,1,,0,0,0
GRID,2,,1,0,0
GRID,3,,1,1,0
GRID,4,,0,1,0
CQUAD4,1,1,1,2,3,4
PSHELL,1,1,0.1,1
MAT1,1,1.0+7,,0.3
SPC1,1,123456,1,2
SPC1,1,12456,3,4
PLOAD2,1,100.,1
FORCE,2,3,0,10.,0.,0.,1.
ENDDATA
"#;
        let m = parse_with_base(deck, None).unwrap();
        assert!(m.independent_steps);
        let out = solve(m).unwrap();
        assert_eq!(out.cases.len(), 2);
        let i3 = out.model.node_index(3).unwrap();
        assert!(
            out.cases[0].u[i3][2] < -1e-6,
            "positive PLOAD2 should push against +Z, uz={}",
            out.cases[0].u[i3][2]
        );
        assert!(out.cases[1].u[i3][2] > 1e-8);
    }

    #[test]
    fn recognized_cards_match_manifest() {
        use crate::mystran_manifest::{names_with, CardStatus};
        let mut got = RECOGNIZED_BULK.to_vec();
        got.sort_unstable();
        let mut expect = names_with(CardStatus::Implemented);
        expect.extend(names_with(CardStatus::MystranBug));
        expect.sort_unstable();
        assert_eq!(got, expect);
    }

    #[test]
    fn declined_card_is_an_error() {
        let deck = "SOL 1\nCEND\nBEGIN BULK\nCUSERIN,1\nGRID,1,,0,0,0\nENDDATA\n";
        let err = parse_with_base(deck, None).unwrap_err();
        assert!(err.to_string().contains("CUSERIN"), "{err}");
    }

    #[test]
    fn crod_matches_checked_in_golden() {
        let deck = r#"
SOL 101
CEND
SPC = 1
LOAD = 1
BEGIN BULK
GRID,1,,0.,0.,0.
GRID,2,,10.,0.,0.
CROD,1,1,1,2
PROD,1,1,2.0
MAT1,1,210000.,,0.3
SPC1,1,123456,1
FORCE,1,2,,1000.,1.,0.,0.
ENDDATA
"#;
        let out = solve(parse_with_base(deck, None).unwrap()).unwrap();
        let golden = include_str!("../tests/mystran/crod_tip.disp");
        for line in golden.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut p = line.split_whitespace();
            let id: i32 = p.next().unwrap().parse().unwrap();
            let exp: [f64; 3] = [
                p.next().unwrap().parse().unwrap(),
                p.next().unwrap().parse().unwrap(),
                p.next().unwrap().parse().unwrap(),
            ];
            let u = out.u[out.model.node_index(id).unwrap()];
            for k in 0..3 {
                let tol = 1e-4 * exp[k].abs() + 1e-8;
                assert!((u[k] - exp[k]).abs() <= tol, "grid {id} u{k}={} exp={}", u[k], exp[k]);
            }
        }
    }

    fn rod_ux(deck: &str) -> f64 {
        let out = solve(parse_with_base(deck, None).unwrap()).unwrap();
        out.u[out.model.node_index(2).unwrap()][0]
    }

    #[test]
    fn syntax_paths_keep_the_rod_exact() {
        let expect = 1000.0 * 10.0 / (210000.0 * 2.0);
        let free = r#"
SOL 101
CEND
ECHO = UNSORT
DISPLACEMENT = ALL
SPC = 1
LOAD = 1
BEGIN BULK
GRID,1,,0.,0.,
,0.
GRID,2,,10.,0.,0.
CROD,1,1,1,2
PROD,1,1,2.0
MAT1,1,2.1+5,,0.3
SPC1,1,123456,1,THRU,1
FORCE,1,2,,1.0+3,1.,0.,0.
ENDDATA
"#;
        assert!((rod_ux(free) - expect).abs() < 1e-8);

        let mut wide = String::from("SOL 101\nCEND\nECHO = NONE\nSPC = 1\nLOAD = 1\nBEGIN BULK\n");
        wide.push_str(&format!(
            "{:<8}{:<16}{:<16}{:<16}{:<16}*\n",
            "GRID*", "1", "0", "0.0", "0.0"
        ));
        wide.push_str(&format!("{:<8}{:<16}\n", "*", "0.0"));
        wide.push_str(&format!(
            "{:<8}{:<16}{:<16}{:<16}{:<16}*\n",
            "GRID*", "2", "0", "10.0", "0.0"
        ));
        wide.push_str(&format!("{:<8}{:<16}\n", "*", "0.0"));
        wide.push_str("CROD,1,1,1,2\nPROD,1,1,2.0\nMAT1,1,210000.,,0.3\n");
        wide.push_str("SPC1,1,123456,1\nFORCE,1,2,,1000.,1.,0.,0.\nENDDATA\n");
        assert!((rod_ux(&wide) - expect).abs() < 1e-8, "wide ux");

        let mut cont = String::from("SOL 1\nCEND\nSPC=1\nLOAD=1\nBEGIN BULK\n");
        let head = format!("{:<8}{:<8}{:<8}{:<8}{:<8}", "GRID", "2", "", "10.0", "0.0");
        cont.push_str(&format!("{head:<72}+G2\n"));
        cont.push_str(&format!("{:<8}{:<8}\n", "+G2", "0.0"));
        cont.push_str("GRID,1,,0.,0.,0.\nCROD,1,1,1,2\nPROD,1,1,2.\nMAT1,1,210000.,,0.3\n");
        cont.push_str("SPC1,1,123456,1\nFORCE,1,2,,1000.,1.,0.,0.\nENDDATA\n");
        assert!((rod_ux(&cont) - expect).abs() < 1e-8, "cont ux {}", rod_ux(&cont));
    }

    #[test]
    fn missing_grid_id_is_an_error() {
        let deck = "SOL 1\nCEND\nBEGIN BULK\nGRID,,,,0,0,0\nENDDATA\n";
        assert!(parse_with_base(deck, None).is_err());
    }

    #[test]
    fn include_is_expanded() {
        let dir = std::env::temp_dir().join(format!("axia-inc-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("mat.bdf"), "MAT1,1,210000.,,0.3\n").unwrap();
        let deck = r#"
SOL 101
CEND
SPC = 1
LOAD = 1
BEGIN BULK
INCLUDE 'mat.bdf'
GRID,1,,0.,0.,0.
GRID,2,,10.,0.,0.
CROD,1,1,1,2
PROD,1,1,2.0
SPC1,1,123456,1
FORCE,1,2,,1000.,1.,0.,0.
ENDDATA
"#;
        let out = solve(parse_with_base(deck, Some(&dir)).unwrap()).unwrap();
        let ux = out.u[out.model.node_index(2).unwrap()][0];
        let expect = 1000.0 * 10.0 / (210000.0 * 2.0);
        assert!((ux - expect).abs() < 1e-8);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
