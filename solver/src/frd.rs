use crate::model::{Model, Procedure};

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
    o.push_str("    1UVERSION           1.9.1\n");
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

    let ne = model.elements.len() as i32;
    o.push_str(&mesh_block_header('3', ne));
    for el in &model.elements {
        let ty = el.kind.frd_type();
        o.push_str(&format!(" -1{:10}{:5}{:5}{:5}\n", el.id, ty, 0, 1));
        o.push_str(" -2");
        for (i, n) in el.nodes.iter().enumerate() {
            o.push_str(&format!("{n:10}"));
            if (i + 1) % 10 == 0 && i + 1 < el.nodes.len() {
                o.push('\n');
                o.push_str(" -2");
            }
        }
        o.push('\n');
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
}
