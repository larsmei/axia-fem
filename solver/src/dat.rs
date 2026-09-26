use crate::model::{Model, Procedure};

pub fn write_dat(
    model: &Model,
    u: &[[f64; 3]],
    stress_gp: &[(i32, usize, [f64; 6])],
    rf: &[[f64; 3]],
) -> String {
    let mut o = String::new();
    if matches!(model.procedure, Procedure::HeatTransfer { .. }) {
        o.push_str(&format!(
            "\n temperatures (NT) for set NALL and time  {:14.7E}\n\n",
            0.0
        ));
        for (k, &id) in model.node_ids.iter().enumerate() {
            o.push_str(&format!("{:10} {:14.6E}\n", id, u[k][0]));
        }
        o.push_str(&format!(
            "\n heat flow (RFL) for set NALL and time  {:14.7E}\n\n",
            0.0
        ));
        for (k, &id) in model.node_ids.iter().enumerate() {
            let r = rf[k][0];
            if r.abs() > 1e-14 {
                o.push_str(&format!("{:10} {:14.6E}\n", id, r));
            }
        }
        o.push('\n');
        return o;
    }
    o.push_str(&format!(
        "\n displacements (vx,vy,vz) for set NALL and time  {:14.7E}\n\n",
        0.0
    ));
    for (k, &id) in model.node_ids.iter().enumerate() {
        o.push_str(&format!(
            "{:10} {:14.6E} {:14.6E} {:14.6E}\n",
            id, u[k][0], u[k][1], u[k][2]
        ));
    }
    o.push_str(&format!(
        "\n forces (fx,fy,fz) for set NALL and time  {:14.7E}\n\n",
        0.0
    ));
    for (k, &id) in model.node_ids.iter().enumerate() {
        let r = rf[k];
        if r[0].abs() + r[1].abs() + r[2].abs() > 1e-14 {
            o.push_str(&format!(
                "{:10} {:14.6E} {:14.6E} {:14.6E}\n",
                id, r[0], r[1], r[2]
            ));
        }
    }
    o.push_str(&format!(
        "\n stresses (elem, integ.pnt.,sxx,syy,szz,sxy,syz,szx) for set EALL and time  {:14.7E}\n\n",
        0.0
    ));
    for (eid, gp, s) in stress_gp {
        o.push_str(&format!(
            "{:10} {:4} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E}\n",
            eid, gp, s[0], s[1], s[2], s[3], s[4], s[5]
        ));
    }
    o.push('\n');
    o
}

/// Surface values (Lamé σθ(a)) live on the nodes. The element mean alone
/// is too coarse for a single quadratic element through the wall.
pub fn append_nodal_stress(o: &mut String, ids: &[i32], stress: &[[f64; 6]]) {
    if stress.is_empty() {
        return;
    }
    o.push_str("\n stresses (node,sxx,syy,szz,sxy,syz,szx) for set NALL\n\n");
    for (k, &id) in ids.iter().enumerate() {
        let s = stress.get(k).copied().unwrap_or([0.0; 6]);
        o.push_str(&format!(
            "{:10} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E}\n",
            id, s[0], s[1], s[2], s[3], s[4], s[5]
        ));
    }
    o.push('\n');
}

pub fn append_frequencies(o: &mut String, freq: &[f64]) {
    if freq.is_empty() {
        return;
    }
    o.push_str("\n EIGENVALUE OUTPUT\n\n MODE  FREQUENCY (CYCLES/TIME)\n\n");
    for (i, f) in freq.iter().enumerate() {
        o.push_str(&format!("{:10} {:14.6E}\n", i + 1, f));
    }
    o.push('\n');
}
