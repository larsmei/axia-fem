mod analysis;
mod beam;
mod dat;
mod elem;
mod error;
mod frd;
mod inp;
mod linalg;
mod model;
mod quadratic;
mod shell;

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
            let mut v = json!({
                "id": e.id,
                "type": e.kind.ccx_name(),
                "nodes": e.nodes,
            });
            if e.kind.is_beam() {
                if let Ok(sec) = model.beam_section_for(e) {
                    v["secA"] = json!(sec.a);
                    v["secB"] = json!(sec.b);
                    v["n1"] = json!(sec.n1);
                    v["area"] = json!(sec.area);
                }
            }
            if e.kind.is_shell() {
                v["th"] = json!(model.thickness_for(e));
            }
            v
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
        "ndofNode": model.ndof_node(),
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
    let ur: Vec<f64> = out.ur.iter().flatten().copied().collect();
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
    v["ur"] = json!(ur);
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

    fn b32_deck(nel: usize, axial: bool) -> String {
        // L=1000, RECT a=10 (n1=Y), b=20 (n2=Z), E=210000, nu=0.3
        let l = 1000.0;
        let nseg = nel; // B32 elements
        let nnode = 2 * nseg + 1;
        let mut s = String::from("*HEADING\nB32 Timoshenko cantilever\n*NODE\n");
        for i in 0..nnode {
            let x = l * i as f64 / (nnode - 1) as f64;
            s.push_str(&format!("{}, {x}, 0, 0\n", i + 1));
        }
        s.push_str("*ELEMENT, TYPE=B32, ELSET=BEAM\n");
        for e in 0..nseg {
            let n1 = 2 * e + 1;
            let mid = n1 + 1;
            let n2 = n1 + 2;
            s.push_str(&format!("{}, {n1}, {n2}, {mid}\n", e + 1));
        }
        s.push_str(
            "*MATERIAL, NAME=STEEL\n*ELASTIC\n210000, 0.3\n\
             *BEAM SECTION, ELSET=BEAM, MATERIAL=STEEL, SECTION=RECT\n\
             10, 20\n0, 1, 0\n*BOUNDARY\n1, 1, 6\n*STEP\n*STATIC\n*CLOAD\n",
        );
        let tip = nnode as i32;
        if axial {
            s.push_str(&format!("{tip}, 1, 21000\n"));
        } else {
            s.push_str(&format!("{tip}, 3, -100\n"));
        }
        s.push_str("*NODE FILE\nU, RF\n*EL FILE\nS\n*END STEP\n");
        s
    }

    #[test]
    fn parse_b32_beam_section() {
        let m = parse_model(&b32_deck(4, false)).unwrap();
        assert_eq!(m.elements.len(), 4);
        assert!(m.has_beams());
        assert_eq!(m.ndof_node(), 6);
        let sec = m.beam_section_for(&m.elements[0]).unwrap();
        assert!((sec.area - 200.0).abs() < 1e-9);
        assert!((sec.i11 - 10.0 * 20.0_f64.powi(3) / 12.0).abs() < 1e-6);
    }

    #[test]
    fn b32_axial_patch() {
        let out = solve_native(&b32_deck(4, true)).unwrap();
        let tip = out.u.last().copied().unwrap();
        // FL/EA = 21000*1000/(210000*200) = 0.5
        assert!(
            (tip[0] - 0.5).abs() < 1e-4,
            "ux_tip={}, expected 0.5",
            tip[0]
        );
    }

    #[test]
    fn b32_cantilever_timoshenko() {
        let out = solve_native(&b32_deck(8, false)).unwrap();
        let tip = out.u.last().copied().unwrap();
        let e = 210000.0;
        let nu = 0.3;
        let g = e / (2.0 * (1.0 + nu));
        let a = 200.0;
        let i11 = 10.0 * 20.0_f64.powi(3) / 12.0;
        let k = 5.0 / 6.0;
        let p = 100.0;
        let l = 1000.0_f64;
        let euler = p * l.powi(3) / (3.0 * e * i11);
        let shear = p * l / (k * g * a);
        let expect = euler + shear;
        let uz = -tip[2];
        assert!(
            (uz - expect).abs() / expect < 0.04,
            "uz={uz}, Timoshenko={expect} (euler={euler} shear={shear})"
        );
        assert!(uz > euler * 0.98, "deflection smaller than Euler — too stiff");
    }

    #[test]
    fn b32_portal_runs() {
        let inp = r#"
*HEADING
Portal
*NODE
1, 0, 0, 0
2, 0, 100, 0
3, 0, 200, 0
4, 200, 200, 0
5, 400, 200, 0
6, 400, 100, 0
7, 400, 0, 0
8, 0, 50, 0
9, 0, 150, 0
10, 100, 200, 0
11, 300, 200, 0
12, 400, 150, 0
13, 400, 50, 0
*ELEMENT, TYPE=B32, ELSET=F
1, 1, 2, 8
2, 2, 3, 9
3, 3, 4, 10
4, 4, 5, 11
5, 5, 6, 12
6, 6, 7, 13
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*BEAM SECTION, ELSET=F, MATERIAL=STEEL, SECTION=RECT
8, 16
0, 0, 1
*BOUNDARY
1, 1, 6
7, 1, 6
*STEP
*STATIC
*CLOAD
5, 1, 50
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        assert!(out.nfree > 0);
        let mut umax = 0.0f64;
        for p in &out.u {
            umax = umax.max((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt());
        }
        assert!(umax > 1e-4 && umax < 50.0, "umax={umax}");
    }

    fn c3d20_patch() -> String {
        r#"
*HEADING
C3D20 uniaxial patch
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
5, 0, 0, 10
6, 10, 0, 10
7, 10, 10, 10
8, 0, 10, 10
9, 5, 0, 0
10, 10, 5, 0
11, 5, 10, 0
12, 0, 5, 0
13, 5, 0, 10
14, 10, 5, 10
15, 5, 10, 10
16, 0, 5, 10
17, 0, 0, 5
18, 10, 0, 5
19, 10, 10, 5
20, 0, 10, 5
*ELEMENT, TYPE=C3D20, ELSET=SOLID
1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=SOLID, MATERIAL=STEEL
*NSET, NSET=FIXED
1, 4, 5, 8, 12, 16, 17, 20
*BOUNDARY
FIXED, 1, 1
1, 2, 3
4, 3, 3
*STEP
*STATIC
*DLOAD
1, P4, 210
*NODE FILE
U, RF
*EL FILE
S
*END STEP
"#
        .into()
    }

    #[test]
    fn patch_test_c3d20() {
        let out = solve_native(&c3d20_patch()).unwrap();
        let mut ux_loaded = Vec::new();
        for (i, &id) in out.model.node_ids.iter().enumerate() {
            if id == 2 || id == 3 || id == 6 || id == 7 {
                ux_loaded.push(out.u[i][0]);
            }
        }
        let mean: f64 = ux_loaded.iter().sum::<f64>() / ux_loaded.len() as f64;
        assert!((mean - 0.01).abs() < 1e-5, "ux={mean}, expected 0.01");
        let mut sxx = 0.0;
        for s in &out.stress {
            sxx += s[0];
        }
        sxx /= out.stress.len() as f64;
        assert!((sxx - 210.0).abs() < 1.0, "sxx={sxx}, expected 210");
    }

    fn cps8_cantilever() -> String {
        let nx = 8usize;
        let ny = 2usize;
        let lx = 100.0;
        let ly = 10.0;
        let mut s = String::from("*HEADING\nCPS8 cantilever\n*NODE\n");
        // Serendipity: no face-centre nodes (odd,odd).
        let mut next = 1i32;
        let mut id_of = vec![vec![0i32; 2 * nx + 1]; 2 * ny + 1];
        for j in 0..=2 * ny {
            for i in 0..=2 * nx {
                if i % 2 == 1 && j % 2 == 1 {
                    continue;
                }
                let x = lx * i as f64 / (2 * nx) as f64;
                let y = ly * j as f64 / (2 * ny) as f64;
                id_of[j][i] = next;
                s.push_str(&format!("{next}, {x}, {y}, 0\n"));
                next += 1;
            }
        }
        s.push_str("*ELEMENT, TYPE=CPS8, ELSET=PLATE\n");
        let mut eid = 1i32;
        for j in 0..ny {
            for i in 0..nx {
                let i0 = 2 * i;
                let j0 = 2 * j;
                let n1 = id_of[j0][i0];
                let n2 = id_of[j0][i0 + 2];
                let n3 = id_of[j0 + 2][i0 + 2];
                let n4 = id_of[j0 + 2][i0];
                let n5 = id_of[j0][i0 + 1];
                let n6 = id_of[j0 + 1][i0 + 2];
                let n7 = id_of[j0 + 2][i0 + 1];
                let n8 = id_of[j0 + 1][i0];
                s.push_str(&format!("{eid}, {n1}, {n2}, {n3}, {n4}, {n5}, {n6}, {n7}, {n8}\n"));
                eid += 1;
            }
        }
        s.push_str(
            "*MATERIAL, NAME=STEEL\n*ELASTIC\n210000, 0.3\n*SOLID SECTION, ELSET=PLATE, MATERIAL=STEEL\n1.0\n",
        );
        s.push_str("*NSET, NSET=FIX\n");
        for j in 0..=2 * ny {
            if id_of[j][0] != 0 {
                s.push_str(&format!("{},\n", id_of[j][0]));
            }
        }
        s.push_str("*BOUNDARY\nFIX, 1, 2\n*STEP\n*STATIC\n*CLOAD\n");
        let tip: Vec<i32> = (0..=2 * ny)
            .filter_map(|j| {
                let id = id_of[j][2 * nx];
                if id == 0 {
                    None
                } else {
                    Some(id)
                }
            })
            .collect();
        let p = 100.0 / tip.len() as f64;
        for id in tip {
            s.push_str(&format!("{id}, 2, {}\n", -p));
        }
        s.push_str("*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n");
        s
    }

    #[test]
    fn patch_test_cps8() {
        let inp = r#"
*HEADING
CPS8 patch
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
5, 5, 0, 0
6, 10, 5, 0
7, 5, 10, 0
8, 0, 5, 0
*ELEMENT, TYPE=CPS8, ELSET=P
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=P, MATERIAL=STEEL
1.0
*BOUNDARY
1, 1, 1
4, 1, 1
8, 1, 1
1, 2, 2
*STEP
*STATIC
*DLOAD
1, P2, -210
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        assert!((u2 - 0.01).abs() < 1e-5, "ux2={u2}");
    }

    #[test]
    fn cantilever_cps8_near_euler() {
        let out = solve_native(&cps8_cantilever()).unwrap();
        let mut umax = 0.0f64;
        for p in &out.u {
            umax = umax.max((-p[1]).max(0.0));
        }
        // Euler ≈ 1.905; quadratic should be closer than linear CPS4
        assert!(
            umax > 1.6 && umax < 2.2,
            "tip |uy|={umax}, expected ~1.9"
        );
    }

    #[test]
    fn parse_quadratic_types() {
        let m = parse_model(&c3d20_patch()).unwrap();
        assert_eq!(m.elements[0].kind.ccx_name(), "C3D20");
        assert_eq!(m.elements[0].nodes.len(), 20);
    }

    fn s4_cantilever(nx: usize, ny: usize) -> String {
        let lx = 100.0;
        let ly = 10.0;
        let mut s = String::from("*HEADING\nS4R Kragplatte\n*NODE\n");
        let id = |i: usize, j: usize| 1 + i + j * (nx + 1);
        for j in 0..=ny {
            for i in 0..=nx {
                let x = lx * i as f64 / nx as f64;
                let y = ly * j as f64 / ny as f64;
                s.push_str(&format!("{}, {x}, {y}, 0\n", id(i, j)));
            }
        }
        s.push_str("*ELEMENT, TYPE=S4R, ELSET=PLATE\n");
        let mut e = 1i32;
        for j in 0..ny {
            for i in 0..nx {
                let n0 = id(i, j);
                let n1 = id(i + 1, j);
                let n2 = id(i + 1, j + 1);
                let n3 = id(i, j + 1);
                s.push_str(&format!("{e}, {n0}, {n1}, {n2}, {n3}\n"));
                e += 1;
            }
        }
        s.push_str(
            "*MATERIAL, NAME=STEEL\n*ELASTIC\n210000, 0.0\n*SHELL SECTION, ELSET=PLATE, MATERIAL=STEEL\n1.0\n",
        );
        s.push_str("*NSET, NSET=FIX\n");
        for j in 0..=ny {
            s.push_str(&format!("{},\n", id(0, j)));
        }
        s.push_str("*BOUNDARY\nFIX, 1, 6\n*STEP\n*STATIC\n*CLOAD\n");
        let p = 1.0 / (ny + 1) as f64;
        for j in 0..=ny {
            s.push_str(&format!("{}, 3, {}\n", id(nx, j), -p));
        }
        s.push_str("*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n");
        s
    }

    #[test]
    fn parse_s4r_shell_section() {
        let m = parse_model(&s4_cantilever(4, 2)).unwrap();
        assert!(m.has_shells());
        assert_eq!(m.ndof_node(), 6);
        assert_eq!(m.elements[0].kind.ccx_name(), "S4R");
        assert!((m.thickness_for(&m.elements[0]) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn s4r_cantilever_near_euler() {
        // L=100, b=10, t=1, E=210000, ν=0, P=1 → δ = PL³/3EI = 1.90476
        let out = solve_native(&s4_cantilever(8, 2)).unwrap();
        let mut umax = 0.0f64;
        for p in &out.u {
            umax = umax.max((-p[2]).max(0.0));
        }
        assert!(
            umax > 1.5 && umax < 2.3,
            "S4R tip |uz|={umax}, expected ~1.905"
        );
    }

    #[test]
    fn s4r_membrane_patch() {
        let inp = r#"
*HEADING
S4 membrane patch
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
*ELEMENT, TYPE=S4, ELSET=P
1, 1, 2, 3, 4
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SHELL SECTION, ELSET=P, MATERIAL=STEEL
1.0
*BOUNDARY
1, 1, 1
4, 1, 1
1, 2, 3
4, 3, 3
2, 3, 3
3, 3, 3
1, 4, 6
2, 4, 6
3, 4, 6
4, 4, 6
*STEP
*STATIC
*CLOAD
2, 1, 1050
3, 1, 1050
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        assert!((u2 - 0.01).abs() < 1e-5, "ux2={u2}");
    }

    #[test]
    fn s4r_ss_plate_pressure() {
        let n = 6usize;
        let mut s = String::from("*HEADING\nSS plate\n*NODE\n");
        let id = |i: usize, j: usize| 1 + i + j * (n + 1);
        for j in 0..=n {
            for i in 0..=n {
                s.push_str(&format!("{}, {}, {}, 0\n", id(i, j), i as f64 * 100.0 / n as f64, j as f64 * 100.0 / n as f64));
            }
        }
        s.push_str("*ELEMENT, TYPE=S4R, ELSET=PLATE\n");
        let mut e = 1;
        for j in 0..n {
            for i in 0..n {
                s.push_str(&format!("{}, {}, {}, {}, {}\n", e, id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)));
                e += 1;
            }
        }
        s.push_str("*MATERIAL, NAME=STEEL\n*ELASTIC\n210000, 0.3\n*SHELL SECTION, ELSET=PLATE, MATERIAL=STEEL\n1.0\n*BOUNDARY\n");
        for i in 0..=n {
            s.push_str(&format!("{}, 3, 3\n", id(i, 0)));
            s.push_str(&format!("{}, 3, 3\n", id(i, n)));
        }
        for j in 1..n {
            s.push_str(&format!("{}, 3, 3\n", id(0, j)));
            s.push_str(&format!("{}, 3, 3\n", id(n, j)));
        }
        s.push_str(&format!("{}, 1, 2\n{}, 2, 2\n", id(0, 0), id(n, 0)));
        s.push_str("*STEP\n*STATIC\n*DLOAD\nEALL, P, -0.01\n*END STEP\n");
        let out = solve_native(&s).unwrap();
        let mut wmax = 0.0f64;
        for p in &out.u {
            wmax = wmax.max((-p[2]).max(0.0));
        }
        // Kirchhoff α q a⁴/D ≈ 0.00406 * 0.01 * 1e8 / 19231 ≈ 0.211
        assert!(
            (wmax - 0.211).abs() < 0.04,
            "wmax={wmax}, expected ~0.211"
        );
    }
}
