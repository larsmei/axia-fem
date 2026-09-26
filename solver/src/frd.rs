use std::collections::HashSet;

use crate::model::{ElemKind, Model, Procedure};

/// Fortran ES12.5 (`1P,E12.5`): sign/space + `d.ddddd` + `E` + `±dd` = 12 chars.
/// Rust `{:.5E}` emits `E-2` (11 chars); CalculiX/Mecway require `E-02`.
fn e12(v: f64) -> String {
    if !v.is_finite() {
        return " 0.00000E+00".to_string();
    }
    if v == 0.0 {
        return " 0.00000E+00".to_string();
    }
    let sign = if v.is_sign_negative() { '-' } else { ' ' };
    let av = v.abs();
    let mut exp = av.log10().floor() as i32;
    let mut mant = av / 10f64.powi(exp);
    if mant >= 10.0 {
        mant /= 10.0;
        exp += 1;
    } else if mant < 1.0 {
        mant *= 10.0;
        exp -= 1;
    }
    let mut scaled = (mant * 1e5).round();
    if scaled >= 1e6 {
        scaled /= 10.0;
        exp += 1;
    }
    let mant_s = format!("{:.5}", scaled / 1e5);
    let exp_s = format!("{exp:+03}");
    let s = format!("{sign}{mant_s}E{exp_s}");
    if s.len() == 12 {
        s
    } else if s.len() > 12 {
        s.chars().take(12).collect()
    } else {
        format!("{s:>12}")
    }
}

fn put_int(buf: &mut [u8], start: usize, width: usize, v: i32) {
    let s = format!("{v:>width$}", width = width);
    let b = s.as_bytes();
    let n = b.len().min(width);
    let off = start + width - n;
    buf[off..off + n].copy_from_slice(&b[b.len() - n..]);
}

/// ccx `fprintf("%5s%1s                  %12d%38d\n", pN, "C", n, 1)` — 74 chars,
/// format flag `1` at column 74 (Mecway `stoi(rec, 74, 75)`).
fn mesh_block_header(kind: char, n: i32) -> String {
    format!("    {kind}C{:18}{n:>12}{:>38}\n", "", 1)
}

/// ccx `frdheader.c` 1PSTEP: 70 chars, counters at cols 25–36 / 37–48 / 49–60.
fn line_1pstep(counter: i32, iinc: i32, istep: i32) -> String {
    let mut t = vec![b' '; 70];
    t[..10].copy_from_slice(b"    1PSTEP");
    put_int(&mut t, 24, 12, counter);
    put_int(&mut t, 36, 12, iinc);
    put_int(&mut t, 48, 12, istep);
    String::from_utf8(t).unwrap() + "\n"
}

/// ccx `frdheader.c` 100CL: 75 chars, ASCII format flag at column 75.
///
/// Field layout (0-based, matching CalculiX 2.22):
/// - [0..7]   `"  100CL"`
/// - [7..12]  `100+iinc` (load-case kode — Mecway keys frames by this)
/// - [12..24] time ES12.5
/// - [24..36] nout
/// - [36..48] description (empty for U/S/RF)
/// - [57]     nmethod flag `'0'` for STATIC
/// - [58..63] iinc
/// - [74]     `'1'` ASCII
///
/// Dataset names live on the following `-4` line, never here. Repeating
/// kode 101/102/103 as dataset types makes Mecway treat every increment
/// as the same load case.
fn line_100cl(iinc: i32, time: f64, nout: i32) -> String {
    let mut t = vec![b' '; 75];
    t[..7].copy_from_slice(b"  100CL");
    put_int(&mut t, 7, 5, 100 + iinc);
    let ts = e12(time);
    t[12..24].copy_from_slice(ts.as_bytes());
    put_int(&mut t, 24, 12, nout);
    t[57] = b'0';
    put_int(&mut t, 58, 5, iinc);
    t[74] = b'1';
    String::from_utf8(t).unwrap() + "\n"
}

fn write_minus2(o: &mut String, nodes: &[i32]) {
    const PER: usize = 10;
    for chunk in nodes.chunks(PER) {
        o.push_str(" -2");
        for n in chunk {
            o.push_str(&format!("{n:10}"));
        }
        o.push('\n');
    }
}

/// ccx `frd.c` reorders midsides for quadratic bricks/wedges so cgx/Mecway
/// see FAM/he20 numbering, not Abaqus INP order.
///
/// C3D20: INP 1–12, 13–16 top, 17–20 vertical → FRD 1–12, 17–20, 13–16.
/// C3D15: INP 1–9, 10–12 top, 13–15 vertical → FRD 1–9, 13, 14, 15, 10–12.
fn frd_nodes(kind: ElemKind, nodes: &[i32]) -> Vec<i32> {
    match kind {
        ElemKind::Hex20 | ElemKind::Hex20R if nodes.len() >= 20 => {
            let n = nodes;
            vec![
                n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[16],
                n[17], n[18], n[19], n[12], n[13], n[14], n[15],
            ]
        }
        ElemKind::Wedge15 if nodes.len() >= 15 => {
            let n = nodes;
            vec![
                n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[12], n[13], n[14], n[9],
                n[10], n[11],
            ]
        }
        // Internal B32 order is Abaqus (end, end, mid). cgx type 12 is end, mid, end.
        ElemKind::Beam32 if nodes.len() >= 3 => vec![nodes[0], nodes[2], nodes[1]],
        _ => nodes.to_vec(),
    }
}

/// Nodes that belong to a mesh element (ccx `inum>0`). Pretension dummy
/// nodes and MASS-only nodes are omitted from 2C and from result blocks.
fn frd_output_nodes(model: &Model) -> Vec<(i32, usize)> {
    let mut used = HashSet::new();
    for el in &model.elements {
        if el.kind.is_point() {
            continue;
        }
        used.extend(el.nodes.iter().copied());
    }
    for p in &model.pretensions {
        used.remove(&p.dummy);
    }
    if used.is_empty() {
        return model
            .node_ids
            .iter()
            .copied()
            .enumerate()
            .map(|(i, id)| (id, i))
            .collect();
    }
    model
        .node_ids
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(i, id)| used.contains(&id).then_some((id, i)))
        .collect()
}

pub struct FrdFrame {
    pub time: f64,
    pub iinc: i32,
    pub u: Vec<[f64; 3]>,
    pub stress: Vec<[f64; 6]>,
    pub rf: Vec<[f64; 3]>,
    pub strain: Vec<[f64; 6]>,
    pub peeq: Vec<f64>,
}

pub fn write_frd(
    model: &Model,
    u: &[[f64; 3]],
    stress: &[[f64; 6]],
    rf: &[[f64; 3]],
    strain: &[[f64; 6]],
    peeq: &[f64],
) -> String {
    write_frd_frames(
        model,
        &[FrdFrame {
            time: 1.0,
            iinc: 1,
            u: u.to_vec(),
            stress: stress.to_vec(),
            rf: rf.to_vec(),
            strain: strain.to_vec(),
            peeq: peeq.to_vec(),
        }],
    )
}

pub fn write_frd_frames(model: &Model, frames: &[FrdFrame]) -> String {
    let mut o = String::with_capacity(1 << 16);
    let heading = if model.heading.is_empty() {
        "Axia"
    } else {
        model.heading.as_str()
    };
    let h: String = heading.chars().take(66).collect();
    o.push_str(&format!("    1C{h}\n"));
    o.push_str("    1UDATE              24.September.2026\n");
    o.push_str("    1UTIME              21:00:00\n");
    o.push_str("    1UHOST              axia\n");
    o.push_str("    1UPGM               Axia FEM\n");
    o.push_str("    1UVERSION           1.31.9\n");
    o.push_str("    1UCODE              CalculiX-compatible Axia FEM\n");

    let out_nodes = frd_output_nodes(model);
    let nn = out_nodes.len() as i32;
    o.push_str(&mesh_block_header('2', nn));
    for &(id, idx) in &out_nodes {
        let c = model.coords[idx];
        o.push_str(&format!(
            " -1{id:10}{}{}{}\n",
            e12(c[0]),
            e12(c[1]),
            e12(c[2])
        ));
    }
    o.push_str(" -3\n");

    // ccx frd.c skips MASS (and other 1-node types). FRD type 11 is a 2-node
    // beam; Mecway reads a second I10 at column 14 and errors
    // "Field of length 10 missing at column 14" on ` -2      8975`.
    let mesh_elems: Vec<_> = model
        .elements
        .iter()
        .filter(|el| !el.kind.is_point())
        .collect();
    let ne = mesh_elems.len() as i32;
    o.push_str(&mesh_block_header('3', ne));
    for el in mesh_elems {
        let ty = el.kind.frd_type();
        o.push_str(&format!(" -1{:10}{:5}{:5}{:5}\n", el.id, ty, 0, 1));
        write_minus2(&mut o, &frd_nodes(el.kind, &el.nodes));
    }
    o.push_str(" -3\n");

    let heat = matches!(model.procedure, Procedure::HeatTransfer { .. });
    let mut icounter = 0i32;
    for (k, fr) in frames.iter().enumerate() {
        let iinc = if fr.iinc > 0 { fr.iinc } else { (k as i32) + 1 };
        write_frame_datasets(&mut o, model, fr, &out_nodes, &mut icounter, iinc, heat);
    }
    o.push_str(" 9999\n");
    o
}

fn write_frame_datasets(
    o: &mut String,
    model: &Model,
    fr: &FrdFrame,
    out_nodes: &[(i32, usize)],
    icounter: &mut i32,
    iinc: i32,
    heat: bool,
) {
    let u = &fr.u;
    let rf = &fr.rf;
    let stress = &fr.stress;
    let strain = &fr.strain;
    let peeq = &fr.peeq;
    let t = if fr.time.is_finite() {
        fr.time
    } else {
        1.0
    };
    // ccx frd.c order for this deck: DISP, STRESS, FORC. 1PSTEP immediately
    // before every 100CL so Mecway sees one load-case header per dataset.
    if heat || model.output_nt {
        let vals: Vec<Vec<f64>> = u.iter().map(|v| vec![v[0]]).collect();
        write_result_block(
            o,
            iinc,
            icounter,
            "NDTEMP",
            1,
            &[("NT", 1, 1, 0)],
            false,
            out_nodes,
            &vals,
            t,
        );
        if model.output_rf {
            let vals: Vec<Vec<f64>> = rf.iter().map(|v| vec![v[0]]).collect();
            write_result_block(
                o,
                iinc,
                icounter,
                "RFL",
                1,
                &[("RFL", 1, 1, 0)],
                false,
                out_nodes,
                &vals,
                t,
            );
        }
    }
    if !heat && model.output_u {
        let vals: Vec<Vec<f64>> = u.iter().map(|v| vec![v[0], v[1], v[2]]).collect();
        write_result_block(
            o,
            iinc,
            icounter,
            "DISP",
            4,
            &[
                ("D1", 2, 1, 0),
                ("D2", 2, 2, 0),
                ("D3", 2, 3, 0),
                ("ALL", 2, 0, 0),
            ],
            true,
            out_nodes,
            &vals,
            t,
        );
    }
    if !heat && model.output_s {
        let vals: Vec<Vec<f64>> = stress
            .iter()
            .map(|v| vec![v[0], v[1], v[2], v[3], v[4], v[5]])
            .collect();
        write_result_block(
            o,
            iinc,
            icounter,
            "STRESS",
            6,
            &[
                ("SXX", 4, 1, 1),
                ("SYY", 4, 2, 2),
                ("SZZ", 4, 3, 3),
                ("SXY", 4, 1, 2),
                ("SYZ", 4, 2, 3),
                ("SZX", 4, 3, 1),
            ],
            false,
            out_nodes,
            &vals,
            t,
        );
    }
    if !heat && model.output_rf {
        let vals: Vec<Vec<f64>> = rf.iter().map(|v| vec![v[0], v[1], v[2]]).collect();
        write_result_block(
            o,
            iinc,
            icounter,
            "FORC",
            4,
            &[
                ("F1", 2, 1, 0),
                ("F2", 2, 2, 0),
                ("F3", 2, 3, 0),
                ("ALL", 2, 0, 0),
            ],
            true,
            out_nodes,
            &vals,
            t,
        );
    }
    if !heat && model.output_e {
        let vals: Vec<Vec<f64>> = strain
            .iter()
            .map(|v| vec![v[0], v[1], v[2], v[3], v[4], v[5]])
            .collect();
        write_result_block(
            o,
            iinc,
            icounter,
            "TOSTRAIN",
            6,
            &[
                ("EXX", 4, 1, 1),
                ("EYY", 4, 2, 2),
                ("EZZ", 4, 3, 3),
                ("EXY", 4, 1, 2),
                ("EYZ", 4, 2, 3),
                ("EZX", 4, 3, 1),
            ],
            false,
            out_nodes,
            &vals,
            t,
        );
    }
    if !heat && peeq.len() == model.node_ids.len() && peeq.iter().any(|v| *v > 0.0) {
        let vals: Vec<Vec<f64>> = peeq.iter().map(|v| vec![*v]).collect();
        write_result_block(
            o,
            iinc,
            icounter,
            "PEEQ",
            1,
            &[("PEEQ", 1, 1, 0)],
            false,
            out_nodes,
            &vals,
            t,
        );
    }
}

fn write_result_block(
    o: &mut String,
    iinc: i32,
    icounter: &mut i32,
    name: &str,
    ncomp_header: i32,
    comps: &[(&str, i32, i32, i32)],
    last_is_all: bool,
    nodes: &[(i32, usize)],
    values: &[Vec<f64>],
    time: f64,
) {
    *icounter += 1;
    let nout = nodes.len() as i32;
    o.push_str(&line_1pstep(*icounter, iinc, 1));
    o.push_str(&line_100cl(iinc, time, nout));
    o.push_str(&format!(" -4  {name:<8}{ncomp_header:5}    1\n"));
    for (i, (cname, typ, a, b)) in comps.iter().enumerate() {
        if last_is_all && i == comps.len() - 1 {
            o.push_str(&format!(
                " -5  {cname:<8}    1{typ:5}{a:5}{b:5}    1ALL\n"
            ));
        } else {
            o.push_str(&format!(" -5  {cname:<8}    1{typ:5}{a:5}{b:5}\n"));
        }
    }
    let nvals = if last_is_all {
        ncomp_header as usize - 1
    } else {
        ncomp_header as usize
    };
    for &(id, idx) in nodes {
        o.push_str(&format!(" -1{id:10}"));
        let empty: Vec<f64> = Vec::new();
        let row = values.get(idx).unwrap_or(&empty);
        for c in 0..nvals {
            let v = row.get(c).copied().unwrap_or(0.0);
            o.push_str(&e12(v));
        }
        o.push('\n');
    }
    o.push_str(" -3\n");
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strip(s: &str) -> &str {
        s.trim_end_matches(['\n', '\r'])
    }

    #[test]
    fn fortran_e12_matches_ccx() {
        assert_eq!(e12(0.06), " 6.00000E-02");
        assert_eq!(e12(0.01), " 1.00000E-02");
        assert_eq!(e12(0.21), " 2.10000E-01");
        assert_eq!(e12(0.0), " 0.00000E+00");
        assert_eq!(e12(1.0), " 1.00000E+00");
        assert_eq!(e12(-1234.5).len(), 12);
        assert_eq!(e12(1.122e7).len(), 12);
        assert!(e12(0.06).chars().nth(8) == Some('E'));
        // two-digit exponent, unlike Rust `{:.5E}` which yields `E-2`
        assert!(e12(0.06).ends_with("E-02"));
        assert_eq!(e12(6e-2).len(), 12);
    }

    #[test]
    fn mesh_headers_are_74_chars_flag_at_col_74() {
        let two = strip(&mesh_block_header('2', 8)).to_string();
        let three = strip(&mesh_block_header('3', 320)).to_string();
        assert_eq!(two.len(), 74, "2C len={}", two.len());
        assert_eq!(three.len(), 74, "3C len={}", three.len());
        assert_eq!(two.as_bytes()[73], b'1');
        assert_eq!(three.as_bytes()[73], b'1');
        assert!(three.starts_with("    3C"));
        // node count occupies the 12-char field after 18 spaces (cols 25–36)
        assert_eq!(&three[24..36], "         320");
        // Mecway's failing line was 62 chars; we must not regress
        assert_ne!(three.len(), 62);
    }

    #[test]
    fn pstep_and_100cl_column_layout() {
        let p = strip(&line_1pstep(1, 1, 1)).to_string();
        assert_eq!(p.len(), 70);
        assert!(p.starts_with("    1PSTEP"));
        assert_eq!(&p[24..36], "           1");
        assert_eq!(&p[36..48], "           1");
        assert_eq!(&p[48..60], "           1");

        let h = strip(&line_100cl(1, 1.0, 8)).to_string();
        assert_eq!(h.len(), 75, "100CL len={}", h.len());
        assert_eq!(h.as_bytes()[74], b'1');
        assert!(h.starts_with("  100CL"));
        assert_eq!(&h[7..12], "  101");
        assert_eq!(&h[58..63], "    1");
        assert_eq!(h.as_bytes()[57], b'0');
        assert!(
            !h.contains("DISP") && !h.contains("STRESS") && !h.contains("FORC"),
            "dataset name belongs on -4, got `{h}`"
        );

        // Byte-identical to CalculiX 2.22 BoltedJoint.frd increment 1 DISP.
        assert_eq!(
            strip(&line_100cl(1, 0.1, 730)),
            "  100CL  101 1.00000E-01         730                     0    1           1"
        );
        assert_eq!(
            strip(&line_100cl(10, 1.0, 730)),
            "  100CL  110 1.00000E+00         730                     0   10           1"
        );
    }

    #[test]
    fn write_frd_mecway_node_and_headers() {
        let inp = r#"
*HEADING
Flacheisen
*NODE
1, 0.06, 0.01, 0.21
2, 1, 0, 0
3, 1, 1, 0
4, 0, 1, 0
5, 0, 0, 1
6, 1, 0, 1
7, 1, 1, 1
8, 0, 1, 1
*ELEMENT, TYPE=C3D8
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=S
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=S
*BOUNDARY
1, 1, 3
*STEP
*STATIC
*CLOAD
2, 1, 1
*END STEP
"#;
        let model = crate::inp::parse(inp).unwrap();
        let z3 = vec![[0.0; 3]; model.node_ids.len()];
        let z6 = vec![[0.0; 6]; model.node_ids.len()];
        let frd = write_frd(&model, &z3, &z6, &z3, &z6, &[]);
        let lines: Vec<&str> = frd.lines().collect();

        let c2 = lines.iter().find(|l| l.starts_with("    2C")).expect("2C");
        let c3 = lines.iter().find(|l| l.starts_with("    3C")).expect("3C");
        assert_eq!(c2.len(), 74, "2C `{c2}`");
        assert_eq!(c3.len(), 74, "3C `{c3}`");
        assert_eq!(c2.as_bytes()[73], b'1');
        assert_eq!(c3.as_bytes()[73], b'1');
        assert_eq!(&c2[24..36], "           8");

        let node = lines
            .iter()
            .find(|l| l.starts_with(" -1") && l.len() >= 13 && l[3..13].trim() == "1")
            .expect("node 1");
        assert!(
            node.starts_with(" -1"),
            "leading space required, got `{node}`"
        );
        assert!(
            node.contains(" 6.00000E-02"),
            "Fortran E-02, got `{node}`"
        );
        assert!(node.contains(" 1.00000E-02"), "{node}");
        assert!(node.contains(" 2.10000E-01"), "{node}");

        assert!(lines.iter().any(|l| *l == " -3"));
        let cl = lines.iter().find(|l| l.starts_with("  100CL")).expect("100CL");
        assert_eq!(cl.len(), 75, "100CL `{cl}`");
        assert_eq!(cl.as_bytes()[74], b'1');
        assert_eq!(&cl[7..12], "  101");
        assert!(!cl.contains("DISP"), "100CL must not carry the dataset name");
    }

    fn dummy_frd(inp: &str) -> String {
        let model = crate::inp::parse(inp).unwrap();
        let z3 = vec![[0.0; 3]; model.node_ids.len()];
        let z6 = vec![[0.0; 6]; model.node_ids.len()];
        write_frd(&model, &z3, &z6, &z3, &z6, &[])
    }

    fn mesh_minus2(frd: &str) -> Vec<String> {
        let mut in_elem = false;
        let mut out = Vec::new();
        for line in frd.lines() {
            if line.starts_with("    3C") {
                in_elem = true;
                continue;
            }
            if in_elem && line == " -3" {
                break;
            }
            if in_elem && line.starts_with(" -2") {
                out.push(line.to_string());
            }
        }
        out
    }

    #[test]
    fn hex20_frd_midsides_match_ccx() {
        // ccx frd.c: kon[0..12], kon[16..20], kon[12..16]
        let n: Vec<i32> = (1..=20).collect();
        assert_eq!(
            frd_nodes(ElemKind::Hex20, &n),
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 17, 18, 19, 20, 13, 14, 15, 16]
        );
        assert_eq!(
            frd_nodes(ElemKind::Hex20R, &n),
            frd_nodes(ElemKind::Hex20, &n)
        );

        let mut nodes = String::new();
        for i in 1..=20 {
            nodes.push_str(&format!("{i}, 0, 0, 0\n"));
        }
        let inp = format!(
            "
*NODE
{nodes}*ELEMENT, TYPE=C3D20
1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
*MATERIAL, NAME=S
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=S
*BOUNDARY
1, 1, 3
*STEP
*STATIC
*CLOAD
2, 1, 1
*END STEP
"
        );
        let minus2 = mesh_minus2(&dummy_frd(&inp));
        assert_eq!(minus2.len(), 2, "{minus2:?}");
        assert_eq!(
            minus2[0],
            " -2         1         2         3         4         5         6         7         8         9        10"
        );
        assert_eq!(
            minus2[1],
            " -2        11        12        17        18        19        20        13        14        15        16"
        );
    }

    #[test]
    fn wedge15_frd_midsides_match_ccx() {
        let n: Vec<i32> = (1..=15).collect();
        assert_eq!(
            frd_nodes(ElemKind::Wedge15, &n),
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 10, 11, 12]
        );

        let mut nodes = String::new();
        for i in 1..=15 {
            nodes.push_str(&format!("{i}, 0, 0, 0\n"));
        }
        let inp = format!(
            "
*NODE
{nodes}*ELEMENT, TYPE=C3D15
1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
*MATERIAL, NAME=S
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=S
*BOUNDARY
1, 1, 3
*STEP
*STATIC
*CLOAD
2, 3, 1
*END STEP
"
        );
        let minus2 = mesh_minus2(&dummy_frd(&inp));
        assert_eq!(minus2.len(), 2, "{minus2:?}");
        assert_eq!(
            minus2[0],
            " -2         1         2         3         4         5         6         7         8         9        13"
        );
        assert_eq!(
            minus2[1],
            " -2        14        15        10        11        12"
        );
    }

    #[test]
    fn mass_omitted_from_frd_mesh_mecway_type11() {
        // Mecway 33 freq-test.inp: S4 then MASS 8723,8975. ccx skips MASS;
        // writing it as FRD type 11 with one node makes Mecway fail:
        // "Field of length 10 missing at column 14 of line  -2      8975"
        let inp = r#"
*NODE
1, 0, 0, 0
2, 1, 0, 0
3, 1, 1, 0
4, 0, 1, 0
8975, -0.291, 0.14, -0.095
*ELEMENT, TYPE=S4
8720, 1, 2, 3, 4
*ELEMENT, TYPE=MASS
8723, 8975
*MATERIAL, NAME=S
*ELASTIC
210000, 0.3
*SHELL SECTION, ELSET=EALL, MATERIAL=S
0.003
*BOUNDARY
1, 1, 6
*STEP
*STATIC
*CLOAD
2, 3, 1
*END STEP
"#;
        let model = crate::inp::parse(inp).unwrap();
        assert_eq!(model.elements.len(), 2);
        assert!(model.elements.iter().any(|e| e.kind.is_point()));
        let frd = dummy_frd(inp);
        let lines: Vec<&str> = frd.lines().collect();

        let c3 = lines
            .iter()
            .find(|l| l.starts_with("    3C"))
            .expect("3C");
        assert_eq!(&c3[24..36], "           1", "3C must count only S4, got `{c3}`");

        let mut in_elem = false;
        let mut elem_ids = Vec::new();
        for line in &lines {
            if line.starts_with("    3C") {
                in_elem = true;
                continue;
            }
            if in_elem && *line == " -3" {
                break;
            }
            if in_elem && line.starts_with(" -1") {
                let id: i32 = line[3..13].trim().parse().expect(line);
                elem_ids.push(id);
            }
            if line.starts_with(" -2") {
                assert!(
                    line.len() >= 23,
                    "Mecway I10 at column 14 missing: `{line}`"
                );
                let nfields = line[3..].split_whitespace().count();
                assert!(nfields >= 2, "too few nodes on `{line}`");
            }
        }
        assert_eq!(elem_ids, vec![8720]);
        assert!(
            !lines.iter().any(|l| *l == " -2      8975"),
            "MASS node must not appear as a 1-node -2 line"
        );
        // MASS-only node 8975 is inum=0 in ccx — omit from 2C.
        let c2 = lines.iter().find(|l| l.starts_with("    2C")).expect("2C");
        assert_eq!(&c2[24..36], "           4", "2C must omit MASS-only node, got `{c2}`");
    }

    #[test]
    fn freq_test_attachment_has_no_short_minus2() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../attachments/freq-test.inp"
        );
        let Ok(inp) = std::fs::read_to_string(path) else {
            return;
        };
        let model = crate::inp::parse(&inp).expect("parse freq-test.inp");
        assert!(
            model.elements.iter().any(|e| e.kind.is_point()),
            "fixture should contain MASS"
        );
        let n = model.node_ids.len();
        let z3 = vec![[0.0; 3]; n];
        let z6 = vec![[0.0; 6]; n];
        let frd = write_frd(&model, &z3, &z6, &z3, &z6, &[]);
        for line in frd.lines() {
            if line.starts_with(" -2") {
                assert!(
                    line.len() >= 23,
                    "Mecway I10 at column 14 missing: `{line}`"
                );
            }
        }
        assert!(
            !frd.lines().any(|l| l == " -2      8975"),
            "exact Mecway failure line still present"
        );
        let c3 = frd.lines().find(|l| l.starts_with("    3C")).unwrap();
        let n_mesh: i32 = c3[24..36].trim().parse().unwrap();
        let n_point = model.elements.iter().filter(|e| e.kind.is_point()).count() as i32;
        assert_eq!(n_mesh, model.elements.len() as i32 - n_point);
    }

    #[test]
    fn write_frd_frames_emits_increment_times() {
        let inp = r#"
*HEADING
inc
*NODE
1, 0, 0, 0
2, 1, 0, 0
3, 1, 1, 0
4, 0, 1, 0
5, 0, 0, 1
6, 1, 0, 1
7, 1, 1, 1
8, 0, 1, 1
*ELEMENT, TYPE=C3D8
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=S
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=S
*BOUNDARY
1, 1, 3
*STEP
*STATIC, DIRECT
0.1, 1
*CLOAD
2, 1, 1
*END STEP
"#;
        let model = crate::inp::parse(inp).unwrap();
        let n = model.node_ids.len();
        let z3 = vec![[0.0; 3]; n];
        let z6 = vec![[0.0; 6]; n];
        let frames: Vec<FrdFrame> = (1..=10)
            .map(|i| FrdFrame {
                time: i as f64 * 0.1,
                iinc: i,
                u: z3.clone(),
                stress: z6.clone(),
                rf: z3.clone(),
                strain: z6.clone(),
                peeq: vec![],
            })
            .collect();
        let frd = write_frd_frames(&model, &frames);
        let lines: Vec<&str> = frd.lines().collect();
        let psteps: Vec<_> = lines
            .iter()
            .copied()
            .filter(|l| l.starts_with("    1PSTEP"))
            .collect();
        // 10 increments × DISP+STRESS+FORC
        assert_eq!(psteps.len(), 30, "1PSTEP before every 100CL");
        let mut n_cl = 0;
        for (i, l) in lines.iter().enumerate() {
            if l.starts_with("  100CL") {
                n_cl += 1;
                assert!(
                    i > 0 && lines[i - 1].starts_with("    1PSTEP"),
                    "100CL must follow 1PSTEP, got prev={}",
                    lines[i.saturating_sub(1)]
                );
                assert!(
                    !l.contains("DISP") && !l.contains("STRESS") && !l.contains("FORC"),
                    "100CL name field must be empty: `{l}`"
                );
            }
        }
        assert_eq!(n_cl, 30);
        let disp: Vec<_> = lines
            .iter()
            .copied()
            .filter(|l| l.starts_with(" -4  DISP"))
            .collect();
        assert_eq!(disp.len(), 10);
        let stress: Vec<_> = lines
            .iter()
            .copied()
            .filter(|l| l.starts_with(" -4  STRESS"))
            .collect();
        assert_eq!(stress.len(), 10);
        let forc: Vec<_> = lines
            .iter()
            .copied()
            .filter(|l| l.starts_with(" -4  FORC"))
            .collect();
        assert_eq!(forc.len(), 10);
        // Dataset order per increment: DISP, STRESS, FORC
        let names: Vec<&str> = lines
            .iter()
            .copied()
            .filter(|l| l.starts_with(" -4  "))
            .take(3)
            .collect();
        assert_eq!(
            names,
            [" -4  DISP        4    1", " -4  STRESS      6    1", " -4  FORC        4    1"]
        );
        let cls: Vec<_> = lines
            .iter()
            .copied()
            .filter(|l| l.starts_with("  100CL"))
            .collect();
        assert!(cls[0].contains(" 1.00000E-01"), "{}", cls[0]);
        assert_eq!(&cls[0][7..12], "  101");
        assert!(cls[27].contains(" 1.00000E+00"), "{}", cls[27]);
        assert_eq!(&cls[27][7..12], "  110");
        assert_eq!(&cls[29][7..12], "  110");
        assert!(frd.contains(" 9999"));
    }
}
