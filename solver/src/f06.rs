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
    s.push_str(&format!(
        "0                                              OUTPUT FOR SUBCASE {n}\n"
    ));
    s.push_str(&format!("0                                              SUBCASE {n}\n"));
    if !case.label.is_empty() {
        s.push_str(&format!("                                               {}\n", case.label));
    }
    if !case.frequencies.is_empty() {
        s.push_str("\n0                                              R E A L   E I G E N V A L U E S\n");
        s.push_str("0     MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES\n");
        for (i, f) in case.frequencies.iter().enumerate() {
            let omega = std::f64::consts::TAU * *f;
            let lambda = omega * omega;
            s.push_str(&format!(
                "     {:6}     {:6}     {}{}{}\n",
                i + 1,
                i + 1,
                es14(lambda),
                es14(omega),
                es14(*f)
            ));
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
    engr_forces(s, model, case);
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
    s.push('\n');
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
    let mut any = false;
    for el in &model.elements {
        if is_rod(el) || is_cbar(model, el) {
            continue;
        }
        any = true;
        break;
    }
    if !any {
        return;
    }
    s.push_str("\n0                                                    S T R E S S E S\n");
    s.push_str("0  ELEMENT  TYPE              SXX           SYY           SZZ           SXY           SYZ           SZX     VON MISES\n");
    for el in &model.elements {
        if is_rod(el) || is_cbar(model, el) {
            continue;
        }
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

fn is_rod(el: &crate::model::Element) -> bool {
    el.kind.is_truss()
}

fn is_cbar(model: &Model, el: &crate::model::Element) -> bool {
    el.kind == crate::model::ElemKind::Beam31
        && model.beam_section_for(el).map(|s| s.cbar).unwrap_or(false)
}

struct BasicDisp {
    u: Vec<[f64; 3]>,
    ur: Vec<[f64; 3]>,
}

fn basic_disp(model: &Model, case: &Subcase) -> BasicDisp {
    let mut u = case.u.clone();
    let mut ur = case.ur.clone();
    if model.output_basic || model.node_transform.is_empty() {
        return BasicDisp { u, ur };
    }
    for (i, id) in model.node_ids.iter().enumerate() {
        let Some(r) = model.node_transform.get(id) else {
            continue;
        };
        if i >= u.len() {
            continue;
        }
        let ul = u[i];
        let rl = ur.get(i).copied().unwrap_or([0.0; 3]);
        for p in 0..3 {
            u[i][p] = r[p][0] * ul[0] + r[p][1] * ul[1] + r[p][2] * ul[2];
            if i < ur.len() {
                ur[i][p] = r[p][0] * rl[0] + r[p][1] * rl[1] + r[p][2] * rl[2];
            }
        }
    }
    BasicDisp { u, ur }
}

fn alpha_dt(model: &Model, el: &crate::model::Element, mat: &crate::model::Material) -> f64 {
    if mat.alpha == 0.0 {
        return 0.0;
    }
    let t = if let Some(&te) = model.elem_temp.get(&el.id) {
        te
    } else {
        let mut s = 0.0;
        let mut n = 0.0;
        for id in &el.nodes {
            s += model.temperature_at(*id);
            n += 1.0;
        }
        if n == 0.0 {
            0.0
        } else {
            s / n
        }
    };
    mat.alpha * (t - mat.tref)
}

fn engr_forces(s: &mut String, model: &Model, case: &Subcase) {
    let disp = basic_disp(model, case);
    let mut rods = Vec::new();
    let mut bars = Vec::new();
    for el in &model.elements {
        if is_rod(el) {
            if let Some(row) = rod_row(model, el, &disp) {
                rods.push((el.id, row));
            }
        } else if is_cbar(model, el) {
            if let Some(row) = bar_row(model, el, &disp) {
                bars.push((el.id, row));
            }
        }
    }
    if !rods.is_empty() {
        write_rod_forces(s, &rods);
        write_rod_stress(s, &rods);
    }
    if !bars.is_empty() {
        write_bar_forces(s, &bars);
        write_bar_stress(s, &bars);
    }
}

struct RodRow {
    force: f64,
    torque: f64,
    stress: f64,
    tau: f64,
}

struct BarRow {
    m1a: f64,
    m2a: f64,
    m1b: f64,
    m2b: f64,
    v1: f64,
    v2: f64,
    axial: f64,
    torque: f64,
    sa: f64,
}

fn rod_row(model: &Model, el: &crate::model::Element, disp: &BasicDisp) -> Option<RodRow> {
    if el.nodes.len() < 2 {
        return None;
    }
    let ia = *model.id_to_index.get(&el.nodes[0])?;
    let ib = *model.id_to_index.get(&el.nodes[1])?;
    let a = *model.coords.get(ia)?;
    let b = *model.coords.get(ib)?;
    let mut d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if len < 1e-18 {
        return None;
    }
    d[0] /= len;
    d[1] /= len;
    d[2] /= len;
    let ua = disp.u.get(ia).copied().unwrap_or([0.0; 3]);
    let ub = disp.u.get(ib).copied().unwrap_or([0.0; 3]);
    let du = (ub[0] - ua[0]) * d[0] + (ub[1] - ua[1]) * d[1] + (ub[2] - ua[2]) * d[2];
    let mat = model.material_for(el).ok()?;
    let area = model.thickness_for(el);
    if area.abs() < 1e-30 || mat.e <= 0.0 {
        return None;
    }
    let nth = mat.e * area * alpha_dt(model, el, &mat);
    let force = mat.e * area * du / len - nth;
    let stress = force / area;
    Some(RodRow {
        force,
        torque: 0.0,
        stress,
        tau: 0.0,
    })
}

fn bar_row(model: &Model, el: &crate::model::Element, disp: &BasicDisp) -> Option<BarRow> {
    if el.nodes.len() < 2 {
        return None;
    }
    let ia = *model.id_to_index.get(&el.nodes[0])?;
    let ib = *model.id_to_index.get(&el.nodes[1])?;
    let xyz = [*model.coords.get(ia)?, *model.coords.get(ib)?];
    let mut ue = [0.0; 12];
    let ua = disp.u.get(ia).copied().unwrap_or([0.0; 3]);
    let ub = disp.u.get(ib).copied().unwrap_or([0.0; 3]);
    let ra = disp.ur.get(ia).copied().unwrap_or([0.0; 3]);
    let rb = disp.ur.get(ib).copied().unwrap_or([0.0; 3]);
    ue[0..3].copy_from_slice(&ua);
    ue[3..6].copy_from_slice(&ra);
    ue[6..9].copy_from_slice(&ub);
    ue[9..12].copy_from_slice(&rb);
    let mat = model.material_for(el).ok()?;
    let sec = model.beam_section_for(el).ok()?;
    let g = crate::beam::cbar_engr_forces(&xyz, &ue, mat.e, mat.nu, &sec, alpha_dt(model, el, &mat)).ok()?;
    let sa = if sec.area.abs() > 1e-30 { g.axial / sec.area } else { 0.0 };
    Some(BarRow {
        m1a: g.m1a,
        m2a: g.m2a,
        m1b: g.m1b,
        m2b: g.m2b,
        v1: g.v1,
        v2: g.v2,
        axial: g.axial,
        torque: g.torque,
        sa,
    })
}

fn write_rod_forces(s: &mut String, rows: &[(i32, RodRow)]) {
    // WRITE_ELEM_ENGR_FORCE formats 301/401/1301/1302/1303, TYPE ROD.
    s.push_str("\n");
    s.push_str(&format!(
        "{:49}E L E M E N T   E N G I N E E R I N G   F O R C E S\n",
        ""
    ));
    s.push_str(&format!(
        "{:53}F O R   E L E M E N T   T Y P E   {:<11}\n",
        "", "ROD"
    ));
    s.push_str("                 Element     Axial        Torque      Element     Axial        Torque      Element     Axial        Torque\n");
    s.push_str("                    ID       Force                       ID       Force                       ID       Force\n");
    let mut max_a = f64::NEG_INFINITY;
    let mut min_a = f64::INFINITY;
    let mut max_t = f64::NEG_INFINITY;
    let mut min_t = f64::INFINITY;
    for chunk in rows.chunks(3) {
        s.push_str("                ");
        for (id, row) in chunk {
            s.push_str(&format!("{:8}{}{}", id, es14(row.force), es14(row.torque)));
            max_a = max_a.max(row.force);
            min_a = min_a.min(row.force);
            max_t = max_t.max(row.torque);
            min_t = min_t.min(row.torque);
        }
        s.push('\n');
    }
    let abs_a = max_a.abs().max(min_a.abs());
    let abs_t = max_t.abs().max(min_t.abs());
    s.push_str("                         ------------- -------------\n");
    s.push_str(&format!("                MAX* :  {}{}\n", es14(max_a), es14(max_t)));
    s.push_str(&format!("                MIN* :  {}{}\n", es14(min_a), es14(min_t)));
    s.push('\n');
    s.push_str(&format!("                ABS* :  {}{}\n", es14(abs_a), es14(abs_t)));
    s.push_str("                *for output set\n");
}

fn write_rod_stress(s: &mut String, rows: &[(i32, RodRow)]) {
    // WRITE_ELEM_STRESSES 301/401/1501 and WRITE_ROD line (margins blank).
    s.push_str("\n");
    s.push_str(&format!(
        "{:29}E L E M E N T   S T R E S S E S   I N   L O C A L   E L E M E N T   C O O R D I N A T E   S Y S T E M\n",
        ""
    ));
    s.push_str(&format!(
        "{:58}F O R   E L E M E N T   T Y P E   {:<11}\n",
        "", "ROD"
    ));
    s.push_str(" Element     Axial       Safety     Torsional     Safety    Element    Axial       Safety     Torsional     Safety\n");
    s.push_str("    ID       Stress      Margin       Stress      Margin       ID      Stress      Margin       Stress      Margin\n");
    let mut max_s = f64::NEG_INFINITY;
    let mut min_s = f64::INFINITY;
    let mut max_t = f64::NEG_INFINITY;
    let mut min_t = f64::INFINITY;
    for chunk in rows.chunks(2) {
        s.push(' ');
        for (id, row) in chunk {
            s.push_str(&format!(
                "{:8}{}{:10} {:1}{}{:10} {:1}",
                id,
                es14(row.stress),
                "",
                "",
                es14(row.tau),
                "",
                ""
            ));
            max_s = max_s.max(row.stress);
            min_s = min_s.min(row.stress);
            max_t = max_t.max(row.tau);
            min_t = min_t.min(row.tau);
        }
        s.push('\n');
    }
    s.push_str("          ------------- ---------  ------------- ---------\n");
    s.push_str(&format!(
        " MAX* :  {}{:10} {}{:10}\n",
        es14(max_s),
        "",
        es14(max_t),
        ""
    ));
    s.push_str(&format!(
        " MIN* :  {}{:10} {}{:10}\n",
        es14(min_s),
        "",
        es14(min_t),
        ""
    ));
    s.push('\n');
    s.push_str(&format!(
        " ABS* :  {}{:10} {}{:10}\n",
        es14(max_s.abs().max(min_s.abs())),
        "",
        es14(max_t.abs().max(min_t.abs())),
        ""
    ));
    s.push_str(" *for output set\n");
}

fn write_bar_forces(s: &mut String, rows: &[(i32, BarRow)]) {
    s.push_str("\n");
    s.push_str(&format!(
        "{:49}E L E M E N T   E N G I N E E R I N G   F O R C E S\n",
        ""
    ));
    s.push_str(&format!(
        "{:61}F O R   E L E M E N T   T Y P E   {:<11}\n",
        "", "BAR"
    ));
    s.push_str("                 Element       Bend-Moment End A           Bend-Moment End B              - Shear -              Axial         Torque\n");
    s.push_str("                    ID       Plane 1       Plane 2       Plane 1       Plane 2      Plane 1       Plane 2        Force\n");
    let mut cols = vec![Vec::new(); 8];
    for (id, row) in rows {
        let v = [
            row.m1a, row.m2a, row.m1b, row.m2b, row.v1, row.v2, row.axial, row.torque,
        ];
        s.push_str("                ");
        s.push_str(&format!("{:8}", id));
        for (k, x) in v.iter().enumerate() {
            s.push_str(&es14(*x));
            cols[k].push(*x);
        }
        s.push('\n');
    }
    s.push_str("                         ------------- ------------- ------------- ------------- ------------- ------------- ------------- -------------\n");
    s.push_str("                MAX* :  ");
    for c in &cols {
        s.push_str(&es14(c.iter().copied().fold(f64::NEG_INFINITY, f64::max)));
    }
    s.push('\n');
    s.push_str("                MIN* :  ");
    for c in &cols {
        s.push_str(&es14(c.iter().copied().fold(f64::INFINITY, f64::min)));
    }
    s.push_str("\n\n");
    s.push_str("                ABS* :  ");
    for c in &cols {
        let mx = c.iter().copied().fold(0.0_f64, |a, v| a.max(v.abs()));
        s.push_str(&es14(mx));
    }
    s.push_str("\n                *for output set\n");
}

fn write_bar_stress(s: &mut String, rows: &[(i32, BarRow)]) {
    s.push_str("\n");
    s.push_str(&format!(
        "{:29}E L E M E N T   S T R E S S E S   I N   L O C A L   E L E M E N T   C O O R D I N A T E   S Y S T E M\n",
        ""
    ));
    s.push_str(&format!(
        "{:58}F O R   E L E M E N T   T Y P E   {:<11}\n",
        "", "BAR"
    ));
    s.push_str(" Element      SA1           SA2           SA3           SA4           Axial        SA-Max        SA-Min      M.S.-T     Torsional\n");
    s.push_str("    ID        SB1           SB2           SB3           SB4          Stress        SB-Max        SB-Min      M.S.-C   Stress/Margin\n");
    for (id, row) in rows {
        // Recovery points C/D/E/F are not stored; SA1–SA4 stay 0, axial is N/A.
        s.push_str(&format!(
            " {:8}{}{}{}{}{}{}{}{:10}{}\n",
            id,
            es14(0.0),
            es14(0.0),
            es14(0.0),
            es14(0.0),
            es14(row.sa),
            es14(row.sa),
            es14(row.sa),
            "",
            es14(0.0),
        ));
        s.push_str(&format!(
            " {:8}{}{}{}{}{:14}{}{}\n",
            "",
            es14(0.0),
            es14(0.0),
            es14(0.0),
            es14(0.0),
            "",
            es14(row.sa),
            es14(row.sa),
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
    // {:.6} can round 9.9999995 up to 10.000000 and break the one-digit mantissa.
    if format!("{mant:.6}").starts_with("10") {
        mant = 1.0;
        exp += 1;
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

    fn near(got: &str, needle: &str) {
        assert!(got.contains(needle), "missing {needle}\n{got}");
    }

    #[test]
    fn crod_force_and_stress_match_the_load() {
        let model = crate::parse_model(
            "SOL 101\nCEND\nSPC=1\nLOAD=1\nBEGIN BULK\n\
             GRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\n\
             CROD,1,1,1,2\nPROD,1,1,2.\nMAT1,1,100.,,0.\n\
             FORCE,1,2,,10.,1.,0.,0.\nSPC1,1,123456,1\nENDDATA\n",
        )
        .unwrap();
        let out = crate::analysis::solve(model.clone()).unwrap();
        let text = write_f06(&model, &out);
        near(&text, "E L E M E N T   E N G I N E E R I N G   F O R C E S");
        near(&text, "T Y P E   ROD");
        near(&text, "1.000000E+01");
        near(&text, "E L E M E N T   S T R E S S E S");
        near(&text, "5.000000E+00");
        assert!((out.u[1][0] - 0.05).abs() < 1e-8, "ux {}", out.u[1][0]);
    }

    #[test]
    fn cbar_tip_forces_are_shear_and_moment() {
        let model = crate::parse_model(
            "SOL 101\nCEND\nSPC=1\nLOAD=1\nBEGIN BULK\n\
             GRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\n\
             CBAR,1,1,1,2,0.,1.,0.\nPBAR,1,1,1.,1.,1.,1.\nMAT1,1,1.,,0.\n\
             FORCE,1,2,,1.,0.,1.,0.\nSPC1,1,123456,1\nENDDATA\n",
        )
        .unwrap();
        let out = crate::analysis::solve(model.clone()).unwrap();
        assert!((out.u[1][1] - 1.0 / 3.0).abs() < 1e-4, "uy {}", out.u[1][1]);
        let text = write_f06(&model, &out);
        near(&text, "T Y P E   BAR");
        near(&text, "Bend-Moment End A");
        near(&text, "1  1.000000E+00  0.000000E+00");
        near(&text, "1.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00");
    }
}
