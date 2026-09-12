use crate::model::{Model, Procedure};

fn e12(v: f64) -> String {
    // Fortran-like 12.5E, always 12 chars. Rust's {:12.5E} is close.
    let s = format!("{:12.5E}", v);
    if s.len() == 12 {
        s
    } else if s.len() > 12 {
        format!("{:11.4E}", v)
    } else {
        format!("{:>12}", s)
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
    o.push_str("    1UVERSION           1.5\n");
    o.push_str("    1UCODE              CalculiX-compatible Axia FEM\n");

    let nn = model.node_ids.len() as i32;
    o.push_str(&format!("    2C{nn:>18}{:>38}\n", 1));
    for (k, &id) in model.node_ids.iter().enumerate() {
        let c = model.coords[k];
        o.push_str(&format!(
            "-1{id:10}{}{}{}\n",
            e12(c[0]),
            e12(c[1]),
            e12(c[2])
        ));
    }
    o.push_str("-3\n");

    let ne = model.elements.len() as i32;
    o.push_str(&format!("    3C{ne:>18}{:>38}\n", 1));
    for el in &model.elements {
        let ty = el.kind.frd_type();
        o.push_str(&format!("-1{:10}{:5}{:5}{:5}\n", el.id, ty, 0, 1));
        o.push_str("-2");
        for (i, n) in el.nodes.iter().enumerate() {
            o.push_str(&format!("{n:10}"));
            if (i + 1) % 10 == 0 && i + 1 < el.nodes.len() {
                o.push('\n');
                o.push_str("-2");
            }
        }
        o.push('\n');
    }
    o.push_str("-3\n");

    o.push_str("    1PSTEP                          1           1           1\n");

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
    // 100CL header matching ccx frdheader roughly:
    // "  100CL" + counter + time + nnodes + itype + name + istep
    let name_pad = format!("{name:<6}");
    o.push_str(&format!(
        "  100CL{kode:5}{:12.5E}{nout:12}{:5}{name_pad:>6}{:5}\n",
        1.0_f64, 0, 1
    ));
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
        o.push_str(&format!("-1{id:10}"));
        let row = &values[k];
        for c in 0..nvals {
            let v = row.get(c).copied().unwrap_or(0.0);
            o.push_str(&e12(v));
        }
        o.push('\n');
    }
    o.push_str("-3\n");
}
