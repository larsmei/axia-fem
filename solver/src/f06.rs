//! MYSTRAN F06 print file.
//!
//! Grid tables follow `WRITE_GRD_PRT_OUTPUTS`: 108-character lines, `ES14.6`
//! fields, headers `D I S P L A C E M E N T S` and `S P C   F O R C E S`.

use crate::analysis::{SolveOutput, Subcase};
use crate::model::Model;

pub fn write_f06(model: &Model, out: &SolveOutput) -> String {
    let mut s = String::new();
    s.push_str("1                                        A X I A   F E M\n");
    s.push_str("                                              MYSTRAN-compatible F06\n");
    s.push_str(&format!("0     {}\n", model.heading));
    s.push('\n');
    if let Some(gid) = model.grdpnt {
        write_grdpnt(&mut s, model, gid);
    }
    if out.cases.is_empty() {
        let one = Subcase {
            label: model.heading.clone(),
            u: out.u.clone(),
            ur: out.ur.clone(),
            rf: out.rf.clone(),
            rm: out.rm.clone(),
            stress: out.stress.clone(),
            von_mises: out.von_mises.clone(),
            frequencies: out.frequencies.clone(),
            buckles: out.buckles.clone(),
            modes: out.modes.clone(),
            mode_ur: out.mode_ur.clone(),
        };
        write_case(&mut s, model, 1, &one);
    } else {
        for (i, case) in out.cases.iter().enumerate() {
            write_case(&mut s, model, i + 1, case);
        }
    }
    s.push_str("0                                              * * * END OF FILE * * *\n");
    s
}

fn write_case(s: &mut String, model: &Model, n: usize, case: &Subcase) {
    s.push_str(&format!("0                                              SUBCASE {n}\n"));
    if !case.label.is_empty() {
        s.push_str(&format!("                                               {}\n", case.label));
    }
    if !case.frequencies.is_empty() {
        s.push_str("\n0                                              R E A L   E I G E N V A L U E S\n");
        s.push_str("0     MODE        FREQUENCY\n");
        s.push_str("                  CYCLES\n");
        for (i, f) in case.frequencies.iter().enumerate() {
            s.push_str(&format!("     {:6}     {}\n", i + 1, es14(*f)));
        }
    }
    if !case.buckles.is_empty() {
        s.push_str("\n0                                          B U C K L I N G   F A C T O R S\n");
        s.push_str("0     MODE        FACTOR\n");
        for (i, f) in case.buckles.iter().enumerate() {
            s.push_str(&format!("     {:6}     {}\n", i + 1, es14(*f)));
        }
    }
    if !case.modes.is_empty() {
        for (i, uu) in case.modes.iter().enumerate() {
            let uur = case.mode_ur.get(i).map(|v| v.as_slice()).unwrap_or(&[]);
            let tag = if !case.buckles.is_empty() {
                format!("MODE {}  FACTOR {}", i + 1, es14(case.buckles.get(i).copied().unwrap_or(0.0)))
            } else {
                format!(
                    "MODE {}  FREQUENCY {}",
                    i + 1,
                    es14(case.frequencies.get(i).copied().unwrap_or(0.0))
                )
            };
            grid_table(
                s,
                &format!("                                                      D I S P L A C E M E N T S"),
                &tag,
                model,
                uu,
                uur,
                false,
            );
        }
    } else {
        grid_table(
            s,
            "                                                      D I S P L A C E M E N T S",
            "(in global coordinate system at each grid)",
            model,
            &case.u,
            &case.ur,
            false,
        );
    }
    grid_table(
        s,
        "                                                         S P C   F O R C E S",
        "(in global coordinate system at each grid)",
        model,
        &case.rf,
        &case.rm,
        true,
    );
    stress_table(s, model, case);
}

fn write_grdpnt(s: &mut String, model: &Model, gid: i32) {
    let origin = if gid == 0 {
        [0.0; 3]
    } else {
        match model.node_index(gid) {
            Ok(i) => model.coords[i],
            Err(_) => {
                s.push_str(&format!("0     GRDPNT: Gitter {gid} fehlt.\n"));
                return;
            }
        }
    };
    let (mass, cg, inertia) = rigid_mass(model, origin);
    s.push_str("\n0                                          R I G I D   B O D Y   M A S S\n");
    s.push_str(&format!("0     GRDPNT {:8}     MASS {}\n", gid, es14(mass)));
    s.push_str(&format!(
        "0     C.G.  {}{}{}\n",
        es14(cg[0]),
        es14(cg[1]),
        es14(cg[2])
    ));
    s.push_str(&format!(
        "0     I(XX,YY,ZZ) {}{}{}\n",
        es14(inertia[0]),
        es14(inertia[1]),
        es14(inertia[2])
    ));
}

fn rigid_mass(model: &Model, origin: [f64; 3]) -> (f64, [f64; 3], [f64; 3]) {
    let mut mass = 0.0;
    let mut mom = [0.0; 3];
    let mut ixx = 0.0;
    let mut iyy = 0.0;
    let mut izz = 0.0;
    let mut add = |m: f64, x: [f64; 3]| {
        if m == 0.0 {
            return;
        }
        mass += m;
        let d = [x[0] - origin[0], x[1] - origin[1], x[2] - origin[2]];
        for k in 0..3 {
            mom[k] += m * d[k];
        }
        ixx += m * (d[1] * d[1] + d[2] * d[2]);
        iyy += m * (d[0] * d[0] + d[2] * d[2]);
        izz += m * (d[0] * d[0] + d[1] * d[1]);
    };
    for el in &model.elements {
        if el.kind != crate::model::ElemKind::Mass {
            continue;
        }
        let Ok(m) = model.mass_for(el) else { continue };
        let Ok(i) = model.node_index(el.nodes[0]) else { continue };
        let arm = model.mass_arms.get(&el.id).copied().unwrap_or([0.0; 3]);
        let x = model.coords[i];
        add(m, [x[0] + arm[0], x[1] + arm[1], x[2] + arm[2]]);
    }
    for link in &model.dof_masses {
        if link.c1 > 2 {
            continue;
        }
        let Ok(i) = model.node_index(link.n1) else { continue };
        add(link.k, model.coords[i]);
    }
    for el in &model.elements {
        if !el.kind.is_truss() {
            continue;
        }
        let Ok(mat) = model.material_for(el) else { continue };
        if mat.density == 0.0 || el.nodes.len() < 2 {
            continue;
        }
        let Ok(ia) = model.node_index(el.nodes[0]) else { continue };
        let Ok(ib) = model.node_index(el.nodes[1]) else { continue };
        let a = model.coords[ia];
        let b = model.coords[ib];
        let len = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2) + (b[2] - a[2]).powi(2)).sqrt();
        let m = mat.density * model.thickness_for(el) * len;
        add(m, [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5, (a[2] + b[2]) * 0.5]);
    }
    let cg = if mass.abs() > 0.0 {
        [origin[0] + mom[0] / mass, origin[1] + mom[1] / mass, origin[2] + mom[2] / mass]
    } else {
        origin
    };
    (mass, cg, [ixx, iyy, izz])
}

fn grid_table(
    s: &mut String,
    title: &str,
    sub: &str,
    model: &Model,
    t: &[[f64; 3]],
    r: &[[f64; 3]],
    totals: bool,
) {
    s.push('\n');
    s.push_str(title);
    s.push('\n');
    s.push_str("1                                             ");
    s.push_str(sub);
    s.push_str("\n\n");
    s.push_str("           GRID     COORD      T1            T2            T3            R1            R2            R3\n");
    s.push_str("                     SYS\n");
    let n = model.node_ids.len().min(t.len());
    let mut acc = [0.0; 6];
    for i in 0..n {
        let row = [
            t[i][0],
            t[i][1],
            t[i][2],
            r.get(i).map(|v| v[0]).unwrap_or(0.0),
            r.get(i).map(|v| v[1]).unwrap_or(0.0),
            r.get(i).map(|v| v[2]).unwrap_or(0.0),
        ];
        for k in 0..6 {
            acc[k] += row[k];
        }
        s.push_str(&grid_row(model.node_ids[i], 0, row));
        s.push('\n');
    }
    if totals {
        s.push_str("                       ------------- ------------- ------------- ------------- ------------- -------------\n");
        s.push_str("    SPC FORCE TOTALS:  ");
        for v in acc {
            s.push_str(&es14(v));
        }
        s.push_str("\n     (for output set)\n");
    }
}

fn stress_table(s: &mut String, model: &Model, case: &Subcase) {
    if case.stress.len() < model.node_ids.len() || model.elements.is_empty() {
        return;
    }
    s.push_str("\n0                                                    S T R E S S E S\n");
    s.push_str("0  ELEMENT  TYPE              SXX           SYY           SZZ           SXY           SYZ           SZX     VON MISES\n");
    for el in &model.elements {
        let mut acc = [0.0; 6];
        let mut n = 0.0;
        let mut vm = 0.0;
        for nid in &el.nodes {
            let Some(&i) = model.id_to_index.get(nid) else {
                continue;
            };
            if i >= case.stress.len() {
                continue;
            }
            for k in 0..6 {
                acc[k] += case.stress[i][k];
            }
            vm += case.von_mises.get(i).copied().unwrap_or(0.0);
            n += 1.0;
        }
        if n == 0.0 {
            continue;
        }
        for k in 0..6 {
            acc[k] /= n;
        }
        vm /= n;
        s.push_str(&format!(
            "  {:8}  {:<8}{}{}{}{}{}{}{}\n",
            el.id,
            el.kind.ccx_name(),
            es14(acc[0]),
            es14(acc[1]),
            es14(acc[2]),
            es14(acc[3]),
            es14(acc[4]),
            es14(acc[5]),
            es14(vm),
        ));
    }
}

fn grid_row(id: i32, cid: i32, v: [f64; 6]) -> String {
    // Fortran 9902: 6X,2(1X,I8),6A14  → 108 columns.
    format!(
        "       {id:8} {cid:8}{}{}{}{}{}{}",
        es14(v[0]),
        es14(v[1]),
        es14(v[2]),
        es14(v[3]),
        es14(v[4]),
        es14(v[5])
    )
}

/// Fortran `ES14.6`: one digit before the decimal, six after, exponent `E+00`.
pub fn es14(v: f64) -> String {
    if !v.is_finite() {
        return format!("{:>14}", "NaN");
    }
    if v == 0.0 {
        return "  0.000000E+00".to_string();
    }
    let neg = v.is_sign_negative();
    let a = v.abs();
    let mut exp = a.log10().floor() as i32;
    let mut mant = a / 10f64.powi(exp);
    if mant >= 9.9999995 {
        mant /= 10.0;
        exp += 1;
    }
    if mant < 1.0 {
        mant *= 10.0;
        exp -= 1;
    }
    let sign = if neg { "-" } else { " " };
    let body = format!("{sign}{mant:.6}E{exp:+03}");
    format!("{body:>14}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn es14_width_and_values() {
        assert_eq!(es14(0.0).len(), 14);
        assert_eq!(es14(1.0).len(), 14);
        assert_eq!(es14(-12.5).len(), 14);
        assert!(es14(1.0).trim().starts_with("1.000000E+00") || es14(1.0).contains("1.000000E+00"));
        assert!(es14(-2.5e-4).contains("2.500000E-04"));
        let row = grid_row(101, 0, [1.0, 0.0, -2.0, 0.0, 0.0, 0.0]);
        assert_eq!(row.len(), 108, "{row} ({})", row.len());
    }
}
