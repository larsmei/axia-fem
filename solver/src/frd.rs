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

/// ccx 100CL: 75 chars, ASCII format flag at C index 74 (column 75).
/// `stoi(rec, 74, 75)` still sees the `1`.
fn line_100cl(kode: i32, time: f64, nout: i32, name: &str) -> String {
    let mut t = vec![b' '; 75];
    t[..7].copy_from_slice(b"  100CL");
    put_int(&mut t, 7, 5, kode);
    let ts = e12(time);
    t[12..24].copy_from_slice(ts.as_bytes());
    put_int(&mut t, 24, 12, nout);
    let nb = name.as_bytes();
    let n = nb.len().min(12);
    t[36..36 + n].copy_from_slice(&nb[..n]);
    put_int(&mut t, 58, 5, kode);
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
        _ => nodes.to_vec(),
    }
}

pub fn write_frd(
    model: &Model,
    u: &[[f64; 3]],
    stress: &[[f64; 6]],
    rf: &[[f64; 3]],
    strain: &[[f64; 6]],
    peeq: &[f64],
) -> String {
    let mut o = String::with_capacity(1 << 16);
    let heading = if model.heading.is_empty() {
        "Axia"
    } else {
        model.heading.as_str()
    };
    let h: String = heading.chars().take(66).collect();
    o.push_str(&format!("    1C{h}\n"));
    o.push_str("    1UDATE              12.September.2026\n");
    o.push_str("    1UTIME              18:00:00\n");
    o.push_str("    1UHOST              axia\n");
    o.push_str("    1UPGM               Axia FEM\n");
    o.push_str("    1UVERSION           1.9.6\n");
    o.push_str("    1UCODE              CalculiX-compatible Axia FEM\n");

    let nn = model.node_ids.len() as i32;
    o.push_str(&mesh_block_header('2', nn));
    for (k, &id) in model.node_ids.iter().enumerate() {
        let c = model.coords[k];
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

    o.push_str(&line_1pstep(1, 1, 1));

    let heat = matches!(model.procedure, Procedure::HeatTransfer { .. });
    if heat || model.output_nt {
        write_result_block(
            &mut o,
            201,
            nn,
            "NDTEMP",
            1,
            &[("NT", 1, 1, 0)],
            false,
            &model.node_ids,
            &u.iter().map(|v| vec![v[0]]).collect::<Vec<_>>(),
        );
        if model.output_rf {
            write_result_block(
                &mut o,
                202,
                nn,
                "RFL",
                1,
                &[("RFL", 1, 1, 0)],
                false,
                &model.node_ids,
                &rf.iter().map(|v| vec![v[0]]).collect::<Vec<_>>(),
            );
        }
    }
    if !heat && model.output_u {
        write_result_block(
            &mut o,
            101,
            nn,
            "DISP",
            4,
            &[
                ("D1", 2, 1, 0),
                ("D2", 2, 2, 0),
                ("D3", 2, 3, 0),
                ("ALL", 2, 0, 0),
            ],
            true,
            &model.node_ids,
            &u.iter().map(|v| vec![v[0], v[1], v[2]]).collect::<Vec<_>>(),
        );
    }
    if !heat && model.output_rf {
        write_result_block(
            &mut o,
            102,
            nn,
            "FORC",
            4,
            &[
                ("F1", 2, 1, 0),
                ("F2", 2, 2, 0),
                ("F3", 2, 3, 0),
                ("ALL", 2, 0, 0),
            ],
            true,
            &model.node_ids,
            &rf.iter()
                .map(|v| vec![v[0], v[1], v[2]])
                .collect::<Vec<_>>(),
        );
    }
    if !heat && model.output_s {
        write_result_block(
            &mut o,
            103,
            nn,
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
            &model.node_ids,
            &stress
                .iter()
                .map(|v| vec![v[0], v[1], v[2], v[3], v[4], v[5]])
                .collect::<Vec<_>>(),
        );
    }
    if !heat && model.output_e {
        write_result_block(
            &mut o,
            104,
            nn,
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
            &model.node_ids,
            &strain
                .iter()
                .map(|v| vec![v[0], v[1], v[2], v[3], v[4], v[5]])
                .collect::<Vec<_>>(),
        );
    }
    if !heat && peeq.len() == model.node_ids.len() && peeq.iter().any(|v| *v > 0.0) {
        write_result_block(
            &mut o,
            108,
            nn,
            "PEEQ",
            1,
            &[("PEEQ", 1, 1, 0)],
            false,
            &model.node_ids,
            &peeq.iter().map(|v| vec![*v]).collect::<Vec<_>>(),
        );
    }
    o.push_str(" 9999\n");
    o
}

fn write_result_block(
    o: &mut String,
    kode: i32,
    nout: i32,
    name: &str,
    ncomp_header: i32,
    comps: &[(&str, i32, i32, i32)],
    last_is_all: bool,
    ids: &[i32],
    values: &[Vec<f64>],
) {
    o.push_str(&line_100cl(kode, 1.0, nout, name));
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
    for (k, &id) in ids.iter().enumerate() {
        o.push_str(&format!(" -1{id:10}"));
        let row = &values[k];
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

        let h = strip(&line_100cl(101, 1.0, 8, "DISP")).to_string();
        assert_eq!(h.len(), 75, "100CL len={}", h.len());
        assert_eq!(h.as_bytes()[74], b'1');
        assert!(h.starts_with("  100CL"));
        assert!(h.contains("DISP"));
    }

    #[test]
    fn write_frd_mecway_node_and_headers() {
        let inp = r#"
*HEADING
Flacheisen
*NODE
519, 0.06, 0.01, 0.21
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

        let node = lines
            .iter()
            .find(|l| l.starts_with(" -1") && l.contains("519"))
            .expect("node 519");
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
}
