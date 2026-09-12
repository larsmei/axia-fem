mod analysis;
mod dat;
mod elem;
mod error;
mod frd;
mod inp;
mod linalg;
mod model;

use serde_json::{json, Value};
use wasm_bindgen::prelude::*;

use crate::error::Result;
use crate::model::Model;

fn mesh_json(model: &Model) -> Value {
    let coords: Vec<f64> = model.coords.iter().flatten().copied().collect();
    let elements: Vec<Value> = model
        .elements
        .iter()
        .map(|e| {
            json!({
                "id": e.id,
                "type": e.kind.ccx_name(),
                "nodes": e.nodes,
            })
        })
        .collect();
    let materials: Vec<Value> = model
        .materials
        .iter()
        .map(|(n, m)| {
            json!({
                "name": n,
                "E": m.e,
                "nu": m.nu,
                "density": m.density,
            })
        })
        .collect();
    json!({
        "heading": model.heading,
        "dim": model.dim,
        "nodeIds": model.node_ids,
        "coords": coords,
        "elements": elements,
        "materials": materials,
        "nnode": model.node_ids.len(),
        "nelem": model.elements.len(),
        "warnings": model.warnings,
    })
}

fn preview_json(inp: &str) -> Result<Value> {
    let model = inp::parse(inp)?;
    let mut v = mesh_json(&model);
    v["ok"] = json!(true);
    v["kind"] = json!("preview");
    Ok(v)
}

fn solve_json(inp: &str) -> Result<Value> {
    let model = inp::parse(inp)?;
    let out = analysis::solve(model)?;
    let mut v = mesh_json(&out.model);
    let u: Vec<f64> = out.u.iter().flatten().copied().collect();
    let s: Vec<f64> = out.stress.iter().flatten().copied().collect();
    let rf: Vec<f64> = out.rf.iter().flatten().copied().collect();
    let mut umax = 0.0f64;
    for p in &out.u {
        let m = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        if m > umax {
            umax = m;
        }
    }
    let mut vmin = f64::MAX;
    let mut vmax = f64::MIN;
    for &m in &out.von_mises {
        vmin = vmin.min(m);
        vmax = vmax.max(m);
    }
    if vmin > vmax {
        vmin = 0.0;
        vmax = 0.0;
    }
    v["ok"] = json!(true);
    v["kind"] = json!("solve");
    v["u"] = json!(u);
    v["stress"] = json!(s);
    v["rf"] = json!(rf);
    v["vonMises"] = json!(out.von_mises);
    v["frd"] = json!(out.frd);
    v["dat"] = json!(out.dat);
    v["stats"] = json!({
        "nnode": out.model.node_ids.len(),
        "nelem": out.model.elements.len(),
        "ndof": out.ndof,
        "nfree": out.nfree,
        "solver": out.solver,
        "iterations": out.iters,
        "residual": out.residual,
        "timeMs": out.time_ms,
        "uMax": umax,
        "vmMin": vmin,
        "vmMax": vmax,
        "nbc": out.model.bcs.len(),
        "ncload": out.model.cloads.len(),
    });
    Ok(v)
}

fn wrap_err(e: impl std::fmt::Display) -> String {
    json!({
        "ok": false,
        "error": e.to_string(),
        "kind": "error",
    })
    .to_string()
}

#[wasm_bindgen]
pub fn preview_inp(inp: &str) -> String {
    match preview_json(inp) {
        Ok(v) => v.to_string(),
        Err(e) => wrap_err(e),
    }
}

#[wasm_bindgen]
pub fn solve_inp(inp: &str) -> String {
    match solve_json(inp) {
        Ok(v) => v.to_string(),
        Err(e) => wrap_err(e),
    }
}

/// Native helpers used by tests.
pub fn parse_model(inp: &str) -> Result<Model> {
    inp::parse(inp)
}

pub fn solve_native(inp: &str) -> Result<analysis::SolveOutput> {
    let model = inp::parse(inp)?;
    analysis::solve(model)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elem::von_mises;

    fn cube_tension() -> String {
        r#"
*HEADING
Uniaxial tension patch test C3D8
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
5, 0, 0, 10
6, 10, 0, 10
7, 10, 10, 10
8, 0, 10, 10
*ELEMENT, TYPE=C3D8, ELSET=SOLID
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=SOLID, MATERIAL=STEEL
*NSET, NSET=FIXED
1, 4, 5, 8
*BOUNDARY
FIXED, 1, 1
1, 2, 3
4, 3, 3
*STEP
*STATIC
*CLOAD
2, 1, 5250
3, 1, 5250
6, 1, 5250
7, 1, 5250
*NODE FILE
U, RF
*EL FILE
S
*END STEP
"#
        .into()
    }

    #[test]
    fn parse_nodes_and_hex() {
        let m = parse_model(&cube_tension()).unwrap();
        assert_eq!(m.node_ids.len(), 8);
        assert_eq!(m.elements.len(), 1);
        assert_eq!(m.dim, 3);
        assert!(m.bcs.len() >= 4);
    }

    #[test]
    fn patch_test_c3d8() {
        let out = solve_native(&cube_tension()).unwrap();
        // ux at x=10 should be FL/EA = 21000*10/(210000*100) = 0.01
        let mut ux_loaded = Vec::new();
        for (i, &id) in out.model.node_ids.iter().enumerate() {
            if id == 2 || id == 3 || id == 6 || id == 7 {
                ux_loaded.push(out.u[i][0]);
            }
        }
        let mean: f64 = ux_loaded.iter().sum::<f64>() / ux_loaded.len() as f64;
        assert!(
            (mean - 0.01).abs() < 1e-6,
            "ux={mean}, expected 0.01"
        );
        let mut sxx = 0.0;
        for s in &out.stress {
            sxx += s[0];
        }
        sxx /= out.stress.len() as f64;
        assert!((sxx - 210.0).abs() < 1e-3, "sxx={sxx}, expected 210");
        assert!(von_mises(&out.stress[0]) > 0.0);
        assert!(out.frd.contains("DISP"));
        assert!(out.frd.contains("STRESS"));
        assert!(out.frd.contains(" 9999"));
        assert!(out.dat.contains("displacements"));
    }

    fn beam_2d() -> String {
        // L=100, h=10, t=1, nx=20, ny=4, E=210000, tip load 100
        let nx = 20usize;
        let ny = 4usize;
        let lx = 100.0;
        let ly = 10.0;
        let mut s = String::from(
            "*HEADING\n2D cantilever CPS4\n*NODE\n",
        );
        let mut nid = 1i32;
        for j in 0..=ny {
            for i in 0..=nx {
                let x = lx * i as f64 / nx as f64;
                let y = ly * j as f64 / ny as f64;
                s.push_str(&format!("{nid}, {x}, {y}, 0\n"));
                nid += 1;
            }
        }
        s.push_str("*ELEMENT, TYPE=CPS4, ELSET=PLATE\n");
        let mut eid = 1i32;
        let col = nx + 1;
        for j in 0..ny {
            for i in 0..nx {
                let n0 = (j * col + i + 1) as i32;
                let n1 = n0 + 1;
                let n2 = n1 + col as i32;
                let n3 = n0 + col as i32;
                s.push_str(&format!("{eid}, {n0}, {n1}, {n2}, {n3}\n"));
                eid += 1;
            }
        }
        s.push_str(
            "*MATERIAL, NAME=STEEL\n*ELASTIC\n210000, 0.3\n*SOLID SECTION, ELSET=PLATE, MATERIAL=STEEL\n1.0\n",
        );
        s.push_str("*NSET, NSET=FIX, GENERATE\n1, ");
        s.push_str(&format!("{}, {}\n", 1 + ny * col, col));
        s.push_str("*NSET, NSET=TIP\n");
        for j in 0..=ny {
            let id = (j * col + nx + 1) as i32;
            s.push_str(&format!("{id},\n"));
        }
        s.push_str("*BOUNDARY\nFIX, 1, 2\n*STEP\n*STATIC\n*CLOAD\n");
        let p = 100.0 / (ny + 1) as f64;
        for j in 0..=ny {
            let id = (j * col + nx + 1) as i32;
            s.push_str(&format!("{id}, 2, {}\n", -p));
        }
        s.push_str("*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n");
        s
    }

    #[test]
    fn cantilever_cps4_order_of_magnitude() {
        let out = solve_native(&beam_2d()).unwrap();
        // Euler: PL^3/(3EI), I=t h^3/12 = 83.333, P=100, L=100, E=210000
        // delta = 100*1e6 / (3*210000*83.333) ≈ 1.905
        let mut umax = 0.0f64;
        for p in &out.u {
            umax = umax.max((-p[1]).max(0.0));
        }
        assert!(
            umax > 1.2 && umax < 2.4,
            "tip |uy|={umax}, expected ~1.9"
        );
    }

    #[test]
    fn nset_generate_and_preview() {
        let inp = cube_tension();
        let j = preview_inp(&inp);
        let v: serde_json::Value = serde_json::from_str(&j).unwrap();
        assert_eq!(v["ok"], true);
        assert_eq!(v["nnode"], 8);
    }
}
