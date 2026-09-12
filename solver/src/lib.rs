mod analysis;
mod axisym;
mod beam;
mod constraint;
mod dat;
mod elem;
mod error;
mod eigen;
mod extra;
mod frd;
mod heat;
mod inp;
mod linalg;
mod material;
mod model;
mod nlgeom;
mod quadratic;
mod shell;

#[cfg(not(target_arch = "wasm32"))]
mod sparse_native;

use serde_json::{json, Value};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

use crate::error::Result;

pub use analysis::SolveOutput;
pub use model::{ElemKind, Model};

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
                "alpha": m.alpha,
                "conductivity": m.conductivity,
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
        "procedure": out.procedure,
        "frequencies": out.frequencies,
        "buckles": out.buckles,
        "nsteps": out.nsteps,
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

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub fn preview_inp(inp: &str) -> String {
    match preview_json(inp) {
        Ok(v) => v.to_string(),
        Err(e) => wrap_err(e),
    }
}

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub fn solve_inp(inp: &str) -> String {
    match solve_json(inp) {
        Ok(v) => v.to_string(),
        Err(e) => wrap_err(e),
    }
}

/// Native helpers used by tests and the CLI.
pub fn parse_model(inp: &str) -> Result<Model> {
    inp::parse(inp)
}

pub fn parse_model_with_base(inp: &str, base: Option<&std::path::Path>) -> Result<Model> {
    inp::parse_with_base(inp, base)
}

pub fn solve_native(inp: &str) -> Result<analysis::SolveOutput> {
    let model = inp::parse(inp)?;
    analysis::solve(model)
}

pub fn solve_native_with_base(
    inp: &str,
    base: Option<&std::path::Path>,
) -> Result<analysis::SolveOutput> {
    let model = inp::parse_with_base(inp, base)?;
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
        #[cfg(not(target_arch = "wasm32"))]
        assert!(
            out.solver.contains("PARDISO") || out.solver.contains("rivrs-sparse"),
            "native backend expected, got {}",
            out.solver
        );
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

    #[test]
    fn patch_test_c3d8i() {
        let inp = cube_tension().replace("C3D8", "C3D8I");
        let out = solve_native(&inp).unwrap();
        let mut ux = Vec::new();
        for (i, &id) in out.model.node_ids.iter().enumerate() {
            if id == 2 || id == 3 || id == 6 || id == 7 {
                ux.push(out.u[i][0]);
            }
        }
        let mean: f64 = ux.iter().sum::<f64>() / ux.len() as f64;
        assert!((mean - 0.01).abs() < 1e-5, "C3D8I ux={mean}");
    }

    #[test]
    fn patch_test_c3d8r() {
        let inp = cube_tension().replace("C3D8", "C3D8R");
        let out = solve_native(&inp).unwrap();
        let mut ux = Vec::new();
        for (i, &id) in out.model.node_ids.iter().enumerate() {
            if id == 2 || id == 3 || id == 6 || id == 7 {
                ux.push(out.u[i][0]);
            }
        }
        let mean: f64 = ux.iter().sum::<f64>() / ux.len() as f64;
        assert!((mean - 0.01).abs() < 2e-4, "C3D8R ux={mean}");
    }

    #[test]
    fn patch_test_c3d6() {
        let inp = r#"
*HEADING
C3D6 uniaxial along prism axis
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 0, 10, 0
4, 0, 0, 10
5, 10, 0, 10
6, 0, 10, 10
*ELEMENT, TYPE=C3D6, ELSET=S
1, 1, 2, 3, 4, 5, 6
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*BOUNDARY
1, 3, 3
2, 3, 3
3, 3, 3
1, 1, 2
2, 2, 2
*STEP
*STATIC
*CLOAD
4, 3, 3500
5, 3, 3500
6, 3, 3500
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u4 = out.u[out.model.node_index(4).unwrap()][2];
        // A = 50, F = 10500, σ = 210, uz = FL/EA = 0.01
        assert!((u4 - 0.01).abs() < 2e-4, "C3D6 uz={u4}");
    }

    #[test]
    fn t3d2_axial() {
        let inp = r#"
*HEADING
T3D2 bar
*NODE
1, 0, 0, 0
2, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*BOUNDARY
1, 1, 3
*STEP
*STATIC
*CLOAD
2, 1, 21000
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        // FL/EA = 21000*1000/(210000*200)=0.5
        assert!((u2 - 0.5).abs() < 1e-6, "T3D2 ux={u2}");
    }

    #[test]
    fn equation_ties_two_nodes() {
        let inp = r#"
*HEADING
two bars + equation
*NODE
1, 0, 0, 0
2, 500, 0, 0
3, 500, 0, 0
4, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
2, 3, 4
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*EQUATION
2
2, 1, 1.0, 3, 1, -1.0
*BOUNDARY
1, 1, 3
3, 2, 3
2, 2, 3
4, 2, 3
*STEP
*STATIC
*CLOAD
4, 1, 21000
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        let u3 = out.u[out.model.node_index(3).unwrap()][0];
        let u4 = out.u[out.model.node_index(4).unwrap()][0];
        assert!((u2 - u3).abs() < 1e-8, "equation u2={u2} u3={u3}");
        // two springs in series, each k=EA/L=210000*200/500=84000, keq=42000
        // u4 = F/keq = 21000/42000 = 0.5
        assert!((u4 - 0.5).abs() < 1e-5, "u4={u4}");
        assert!((u2 - 0.25).abs() < 1e-5, "u2={u2}");
    }

    #[test]
    fn springa_axial() {
        let inp = r#"
*HEADING
spring
*NODE
1, 0, 0, 0
2, 1, 0, 0
*ELEMENT, TYPE=SPRINGA, ELSET=S
1, 1, 2
*SPRING, ELSET=S
1000
*BOUNDARY
1, 1, 3
2, 2, 3
*STEP
*STATIC
*CLOAD
2, 1, 10
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        assert!((u2 - 0.01).abs() < 1e-8, "spring ux={u2}");
    }

    #[test]
    fn include_expands() {
        let dir = std::env::temp_dir().join("axia_include_test");
        let _ = std::fs::create_dir_all(&dir);
        std::fs::write(
            dir.join("mat.inp"),
            "*MATERIAL, NAME=STEEL\n*ELASTIC\n210000, 0.3\n",
        )
        .unwrap();
        let main = format!(
            "*HEADING\ninc\n*NODE\n1,0,0,0\n2,1000,0,0\n*ELEMENT, TYPE=T3D2, ELSET=T\n1,1,2\n*INCLUDE, INPUT=mat.inp\n*SOLID SECTION, ELSET=T, MATERIAL=STEEL\n200\n*BOUNDARY\n1,1,3\n*STEP\n*STATIC\n*CLOAD\n2,1,21000\n*END STEP\n"
        );
        std::fs::write(dir.join("job.inp"), &main).unwrap();
        let text = std::fs::read_to_string(dir.join("job.inp")).unwrap();
        let model = crate::inp::parse_with_base(&text, Some(&dir)).unwrap();
        assert!(model.materials.contains_key("STEEL"));
        let out = crate::analysis::solve(model).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        assert!((u2 - 0.5).abs() < 1e-6);
    }

    #[test]
    fn tie_two_truss_nodes() {
        let inp = r#"
*HEADING
tie
*NODE
1, 0, 0, 0
2, 500, 0, 0
3, 500, 0, 0
4, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
2, 3, 4
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*SURFACE, NAME=SL, TYPE=NODE
2
*SURFACE, NAME=MA, TYPE=NODE
3
*TIE
SL, MA
*BOUNDARY
1, 1, 3
2, 2, 3
3, 2, 3
4, 2, 3
*STEP
*STATIC
*CLOAD
4, 1, 21000
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        let u4 = out.u[out.model.node_index(4).unwrap()][0];
        assert!((u2 - 0.25).abs() < 1e-5, "tie u2={u2}");
        assert!((u4 - 0.5).abs() < 1e-5, "tie u4={u4}");
    }

    #[test]
    fn rigid_body_fixed_ref() {
        let inp = r#"
*HEADING
rigid
*NODE
1, 0, 0, 0
2, 100, 0, 0
3, 0, 100, 0
*NSET, NSET=SLAVES
2, 3
*RIGID BODY, NSET=SLAVES, REF NODE=1
*ELEMENT, TYPE=T3D2, ELSET=T
1, 2, 3
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
10
*BOUNDARY
1, 1, 6
*STEP
*STATIC
*CLOAD
2, 2, 50
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()];
        let u3 = out.u[out.model.node_index(3).unwrap()];
        assert!(u2.iter().all(|v| v.abs() < 1e-10), "u2={u2:?}");
        assert!(u3.iter().all(|v| v.abs() < 1e-10), "u3={u3:?}");
        assert_eq!(out.model.ndof_node(), 6);
    }

    #[test]
    fn transform_local_spc() {
        // Local x = global y. Fix local-x (global y) at node 1, load local-x at node 2.
        let inp = r#"
*HEADING
transform
*NODE
1, 0, 0, 0
2, 0, 1000, 0
*NSET, NSET=ALLN
1, 2
*TRANSFORM, NSET=ALLN, TYPE=R
0, 1, 0, -1, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*BOUNDARY
1, 1, 3
*STEP
*STATIC
*CLOAD
2, 1, 21000
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        // local x = global y, bar along global y → axial in local x
        let u2 = out.u[out.model.node_index(2).unwrap()];
        assert!(
            (u2[1] - 0.5).abs() < 1e-5,
            "transform uy (global)={:?} expected 0.5",
            u2
        );
    }

    #[test]
    fn distributing_coupling_translates() {
        let inp = r#"
*HEADING
rbe3
*NODE
1, 0, 0, 0
2, 100, 0, 0
3, 0, 0, 0
*SURFACE, NAME=S, TYPE=NODE
1, 2
*COUPLING, REF NODE=3, SURFACE=S
*DISTRIBUTING
1, 3
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*BOUNDARY
1, 1, 3
*STEP
*STATIC
*CLOAD
3, 1, 21000
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        // ref u = average of 1 and 2; load on ref goes to surface (RBE3)
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        let u3 = out.u[out.model.node_index(3).unwrap()][0];
        assert!(u2 > 0.0, "u2={u2}");
        assert!((u3 - 0.5 * u2).abs() < 1e-6, "u3={u3} should be avg of 0 and u2={u2}");
    }

    #[test]
    fn frequency_t3d2_axial() {
        let inp = r#"
*HEADING
freq
*NODE
1, 0, 0, 0
2, 1, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
100, 0.0
*DENSITY
1.0
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
1.0
*BOUNDARY
1, 1, 3
2, 2, 3
*STEP
*FREQUENCY
1
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        assert!(!out.frequencies.is_empty(), "no frequencies");
        let f = out.frequencies[0];
        // k=EA/L=100, m=rhoAL/2=0.5, f=sqrt(k/m)/(2π)=2.2508
        assert!(
            (f - 2.2508).abs() / 2.2508 < 0.08,
            "f={f}, expected ~2.251 Hz"
        );
    }

    #[test]
    fn buckle_b31_cantilever() {
        // L=1000, RECT 10×20, Imin=1666.67, E=210000, P=1
        // Euler cantilever π²EI/(4L²) ≈ 863.7
        let nseg = 8usize;
        let l = 1000.0;
        let mut s = String::from("*HEADING\nbuckle\n*NODE\n");
        for i in 0..=nseg {
            s.push_str(&format!("{}, {}, 0, 0\n", i + 1, l * i as f64 / nseg as f64));
        }
        s.push_str("*ELEMENT, TYPE=B31, ELSET=B\n");
        for e in 0..nseg {
            s.push_str(&format!("{}, {}, {}\n", e + 1, e + 1, e + 2));
        }
        s.push_str(
            "*MATERIAL, NAME=STEEL\n*ELASTIC\n210000, 0.3\n\
             *BEAM SECTION, ELSET=B, MATERIAL=STEEL, SECTION=RECT\n\
             10, 20\n0, 0, 1\n*BOUNDARY\n1, 1, 6\n*STEP\n*BUCKLE\n1\n*CLOAD\n",
        );
        s.push_str(&format!("{}, 1, -1\n*END STEP\n", nseg + 1));
        let out = solve_native(&s).unwrap();
        assert!(!out.buckles.is_empty(), "no buckle factors");
        let lam = out.buckles[0].abs();
        let expect = std::f64::consts::PI.powi(2) * 210000.0 * (20.0 * 10.0_f64.powi(3) / 12.0)
            / (4.0 * l.powi(2));
        assert!(
            (lam - expect).abs() / expect < 0.25,
            "λ={lam}, Euler cantilever={expect}"
        );
    }

    #[test]
    fn thermal_t3d2_reaction() {
        let inp = r#"
*HEADING
thermal bar
*NODE
1, 0, 0, 0
2, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*EXPANSION
1e-5
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*BOUNDARY
1, 1, 3
2, 1, 3
*TEMPERATURE
1, 100
2, 100
*STEP
*STATIC
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let r1 = out.rf[out.model.node_index(1).unwrap()][0];
        // N = EA α ΔT = 210000*200*1e-5*100 = 42000, compression → node1 reaction +42000 (against expansion)
        assert!(
            (r1 - 42000.0).abs() / 42000.0 < 0.02,
            "R1x={r1}, expected 42000"
        );
    }

    #[test]
    fn nlgeom_truss_matches_small_strain() {
        let inp = r#"
*HEADING
nlgeom small
*NODE
1, 0, 0, 0
2, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*BOUNDARY
1, 1, 3
2, 2, 3
*STEP, NLGEOM
*STATIC
*CLOAD
2, 1, 21000
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        assert!((u2 - 0.5).abs() < 1e-3, "NLGEOM ux={u2}");
    }

    #[test]
    fn plastic_truss_hardening() {
        let inp = r#"
*HEADING
plastic bar
*NODE
1, 0, 0, 0
2, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*PLASTIC
210, 0.0
420, 0.01
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*BOUNDARY
1, 1, 3
2, 2, 3
*STEP
*STATIC
*CLOAD
2, 1, 50000
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        // σ=250, H=21000, peeq=(250-210)/H=0.0019048, ε=250/E+peeq=0.003095, u=3.095
        assert!(
            (u2 - 3.095).abs() / 3.095 < 0.05,
            "plastic ux={u2}, expected ~3.095"
        );
        assert!(u2 > 1.2, "should exceed elastic u=1.19");
    }

    #[test]
    fn plastic_below_yield_is_elastic() {
        let inp = r#"
*HEADING
elastic plastic
*NODE
1, 0, 0, 0
2, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*PLASTIC
210, 0.0
420, 0.01
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*BOUNDARY
1, 1, 3
2, 2, 3
*STEP
*STATIC
*CLOAD
2, 1, 40000
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        // σ=200 < 210 → u = FL/EA = 0.95238
        assert!((u2 - 0.95238).abs() < 1e-3, "elastic plastic ux={u2}");
        assert!(out.iters >= 1);
    }

    #[test]
    fn nlgeom_mixed_mesh_is_rejected() {
        let inp = r#"
*HEADING
nlgeom mixed
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
*ELEMENT, TYPE=S4, ELSET=S
1, 1, 2, 3, 4
*ELEMENT, TYPE=T3D2, ELSET=T
2, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SHELL SECTION, ELSET=S, MATERIAL=STEEL
1.0
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
1
*BOUNDARY
1, 1, 6
*STEP, NLGEOM
*STATIC
*CLOAD
2, 3, 1
*END STEP
"#;
        let e = match solve_native(inp) {
            Ok(_) => panic!("expected mixed NLGEOM to fail"),
            Err(e) => e.to_string(),
        };
        assert!(
            e.contains("NLGEOM") || e.contains("T3D2") || e.contains("Kontinuum"),
            "unexpected error: {e}"
        );
    }

    #[test]
    fn heat_t3d2_linear() {
        let inp = r#"
*HEADING
heat bar
*NODE
1, 0, 0, 0
2, 500, 0, 0
3, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
2, 2, 3
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*CONDUCTIVITY
50
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
1
*BOUNDARY
1, 11, 11, 0
3, 11, 11, 100
*STEP
*HEAT TRANSFER, STEADY STATE
*NODE FILE
NT
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        assert!(out.procedure.contains("HEAT"));
        let t2 = out.u[out.model.node_index(2).unwrap()][0];
        assert!((t2 - 50.0).abs() < 1e-6, "T2={t2}, expected 50");
        assert!(out.frd.contains("NDTEMP") || out.frd.contains("NT"));
        assert!(out.dat.contains("temperatures"));
    }

    #[test]
    fn heat_c3d8_patch() {
        let inp = r#"
*HEADING
heat hex
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
5, 0, 0, 10
6, 10, 0, 10
7, 10, 10, 10
8, 0, 10, 10
*ELEMENT, TYPE=C3D8, ELSET=S
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*CONDUCTIVITY
25
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*NSET, NSET=COLD
1, 4, 5, 8
*NSET, NSET=HOT
2, 3, 6, 7
*BOUNDARY
COLD, 11, 11, 0
HOT, 11, 11, 100
*STEP
*HEAT TRANSFER, STEADY STATE
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        for (i, &id) in out.model.node_ids.iter().enumerate() {
            let x = out.model.coords[i][0];
            let t = out.u[i][0];
            let expect = 10.0 * x;
            assert!(
                (t - expect).abs() < 1e-6,
                "node {id} T={t}, expected {expect}"
            );
        }
    }

    #[test]
    fn dynamic_sdof_half_period() {
        // k = EA/L = 100, m_lumped free = ρAL/2 = 1 → ω = 10 rad/s
        // u(0)=0.01, v(0)=0 → u(π/10) = -0.01
        let inp = r#"
*HEADING
sdof
*NODE
1, 0, 0, 0
2, 1, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
100, 0.0
*DENSITY
2.0
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
1.0
*BOUNDARY
1, 1, 3
2, 2, 3
*INITIAL CONDITIONS, TYPE=DISPLACEMENT
2, 1, 0.01
*STEP
*DYNAMIC
0.005, 0.3141592653589793
*NODE FILE
U
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        assert!(out.procedure.contains("DYNAMIC"));
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        assert!(
            (u2 + 0.01).abs() < 5e-4,
            "u2={u2}, expected -0.01 at T/2"
        );
        assert!(out.solver.contains("Newmark"), "solver={}", out.solver);
    }

    #[test]
    fn amplitude_ramps_static_end() {
        // At t=0 amplitude is unused in STATIC; load is the CLOAD magnitude.
        let inp = r#"
*HEADING
amp static
*NODE
1, 0, 0, 0
2, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*BOUNDARY
1, 1, 3
2, 2, 3
*AMPLITUDE, NAME=RAMP
0, 0
1, 1
*STEP
*STATIC
*CLOAD, AMPLITUDE=RAMP
2, 1, 21000
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        // STATIC uses amp(t=0)=0 → u=0
        assert!(u2.abs() < 1e-8, "static amp(0) ux={u2}");
    }

    #[test]
    fn patch_test_c3d15() {
        // confined uniaxial strain εz=0.001 → σz = E(1-ν)/((1+ν)(1-2ν))*εz
        let inp = r#"
*HEADING
C3D15 confined patch
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 0, 10, 0
4, 0, 0, 10
5, 10, 0, 10
6, 0, 10, 10
7, 5, 0, 0
8, 5, 5, 0
9, 0, 5, 0
10, 5, 0, 10
11, 5, 5, 10
12, 0, 5, 10
13, 0, 0, 5
14, 10, 0, 5
15, 0, 10, 5
*ELEMENT, TYPE=C3D15, ELSET=S
1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*BOUNDARY
1, 1, 2
2, 1, 2
3, 1, 2
4, 1, 2
5, 1, 2
6, 1, 2
7, 1, 2
8, 1, 2
9, 1, 2
10, 1, 2
11, 1, 2
12, 1, 2
13, 1, 2
14, 1, 2
15, 1, 2
1, 3, 3, 0.0
2, 3, 3, 0.0
3, 3, 3, 0.0
7, 3, 3, 0.0
8, 3, 3, 0.0
9, 3, 3, 0.0
13, 3, 3, 0.005
14, 3, 3, 0.005
15, 3, 3, 0.005
4, 3, 3, 0.01
5, 3, 3, 0.01
6, 3, 3, 0.01
10, 3, 3, 0.01
11, 3, 3, 0.01
12, 3, 3, 0.01
*STEP
*STATIC
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let sz: f64 = out.stress.iter().map(|s| s[2]).sum::<f64>() / out.stress.len() as f64;
        let expect = 210000.0 * 0.7 / (1.3 * 0.4) * 0.001;
        assert!(
            (sz - expect).abs() / expect < 0.05,
            "C3D15 σz={sz}, expected {expect}"
        );
    }

    #[test]
    fn patch_test_cax4_lame() {
        // infinite cylinder, internal pressure. 4 CAX4 through the wall.
        let mut inp = String::from(
            "*HEADING\nCAX4 Lame\n*NODE\n",
        );
        let a = 10.0;
        let b = 20.0;
        let nr = 4;
        for i in 0..=nr {
            let r = a + (b - a) * i as f64 / nr as f64;
            let n0 = 1 + 2 * i;
            let n1 = n0 + 1;
            inp.push_str(&format!("{n0}, {r}, 0\n{n1}, {r}, 2\n"));
        }
        inp.push_str("*ELEMENT, TYPE=CAX4, ELSET=S\n");
        for i in 0..nr {
            let n0 = 1 + 2 * i;
            let n1 = n0 + 2;
            let n2 = n1 + 1;
            let n3 = n0 + 1;
            inp.push_str(&format!("{}, {}, {}, {}, {}\n", i + 1, n0, n1, n2, n3));
        }
        inp.push_str(
            "*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*BOUNDARY
",
        );
        for i in 0..=nr {
            let n0 = 1 + 2 * i;
            let n1 = n0 + 1;
            inp.push_str(&format!("{n0}, 2, 2\n{n1}, 2, 2\n"));
        }
        inp.push_str(
            "*STEP
*STATIC
*DLOAD
1, P4, 10
*END STEP
",
        );
        // face 4 of first element: nodes n3-n0 = inner edge (3-0 of elem 1 = 2,1) wait
        // elem 1 nodes: n0=1 (r=a,z=0), n1=3 (r=a+dr,z=0), n2=4 (r=a+dr,z=2), n3=2 (r=a,z=2)
        // edges: P1=1-2 (bottom), P2=2-3 (outer), P3=3-4 (top), P4=4-1 (inner) → P4 is inner. Good.
        let out = solve_native(&inp).unwrap();
        let ur = out.u[out.model.node_index(1).unwrap()][0];
        let pin = 10.0;
        let aa = a * a;
        let bb = b * b;
        let a_const = pin * aa / (bb - aa);
        let b_const = pin * aa * bb / (bb - aa);
        let s_th = a_const + b_const / aa;
        let s_r = a_const - b_const / aa;
        let s_z = 2.0 * 0.3 * a_const;
        let eth = (s_th - 0.3 * s_r - 0.3 * s_z) / 210000.0;
        let expect = a * eth;
        assert!(
            (ur - expect).abs() / expect.abs() < 0.12,
            "CAX4 ur={ur}, Lame {expect}"
        );
    }

    #[test]
    fn patch_test_m3d4() {
        let inp = r#"
*HEADING
M3D4 membrane patch
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
*ELEMENT, TYPE=M3D4, ELSET=M
1, 1, 2, 3, 4
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.0
*MEMBRANE SECTION, ELSET=M, MATERIAL=STEEL
1.0
*BOUNDARY
1, 1, 1
4, 1, 1
1, 2, 2
2, 2, 2
1, 3, 3
2, 3, 3
3, 3, 3
4, 3, 3
2, 1, 1, 0.01
3, 1, 1, 0.01
*STEP
*STATIC
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        assert!((u2 - 0.01).abs() < 1e-8, "M3D4 ux={u2}");
        let sxx = out.stress[out.model.node_index(3).unwrap()][0];
        assert!((sxx - 210.0).abs() < 1.0, "M3D4 sxx={sxx}");
    }

    #[test]
    fn heat_c3d20_linear() {
        let inp = r#"
*HEADING
heat hex20
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
*ELEMENT, TYPE=C3D20, ELSET=S
1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*CONDUCTIVITY
25
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*BOUNDARY
1, 11, 11, 0
4, 11, 11, 0
5, 11, 11, 0
8, 11, 11, 0
12, 11, 11, 0
16, 11, 11, 0
17, 11, 11, 0
20, 11, 11, 0
2, 11, 11, 100
3, 11, 11, 100
6, 11, 11, 100
7, 11, 11, 100
10, 11, 11, 100
14, 11, 11, 100
18, 11, 11, 100
19, 11, 11, 100
*STEP
*HEAT TRANSFER, STEADY STATE
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let t9 = out.u[out.model.node_index(9).unwrap()][0];
        assert!((t9 - 50.0).abs() < 1e-4, "C3D20 mid-edge T={t9}");
    }

    #[test]
    fn buckle_c3d8_column_positive() {
        let inp = r#"
*HEADING
C3D8 column buckle
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
5, 0, 0, 100
6, 10, 0, 100
7, 10, 10, 100
8, 0, 10, 100
*ELEMENT, TYPE=C3D8, ELSET=S
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*BOUNDARY
1, 1, 3
2, 2, 3
3, 3, 3
4, 1, 3
2, 1, 1
4, 2, 2
*STEP
*BUCKLE
1
*CLOAD
5, 3, -250
6, 3, -250
7, 3, -250
8, 3, -250
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        assert!(!out.buckles.is_empty(), "no buckle factors");
        assert!(
            out.buckles[0].is_finite() && out.buckles[0] > 0.0,
            "lambda={}",
            out.buckles[0]
        );
    }

    #[test]
    fn two_static_steps_accumulate_cload() {
        let inp = r#"
*HEADING
two steps
*NODE
1, 0, 0, 0
2, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*BOUNDARY
1, 1, 3
*STEP
*STATIC
*CLOAD
2, 1, 10500
*END STEP
*STEP
*STATIC
*CLOAD
2, 1, 10500
*END STEP
"#;
        let model = parse_model(inp).unwrap();
        assert_eq!(model.steps.len(), 2, "steps={}", model.steps.len());
        let out = solve_native(inp).unwrap();
        assert_eq!(out.nsteps, 2);
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        assert!((u2 - 0.5).abs() < 1e-6, "two-step ux={u2}");
    }

    #[test]
    fn controls_sets_maxiter() {
        let inp = r#"
*HEADING
controls
*NODE
1, 0, 0, 0
2, 1, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
100, 0
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
1
*BOUNDARY
1, 1, 3
*CONTROLS, MAXITER=7, RTOL=1e-10
*STEP
*STATIC
*CLOAD
2, 1, 1
*END STEP
"#;
        let model = parse_model(inp).unwrap();
        assert_eq!(model.max_newton, 7);
        assert!((model.newton_tol - 1e-10).abs() < 1e-20);
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        assert!((u2 - 0.01).abs() < 1e-8, "ux={u2}");
    }

    #[test]
    fn nlgeom_c3d8_small_strain_patch() {
        let mut inp = cube_tension();
        inp = inp.replace("*STEP\n*STATIC", "*STEP, NLGEOM\n*STATIC");
        let out = solve_native(&inp).unwrap();
        let mut ux_loaded = Vec::new();
        for (i, &id) in out.model.node_ids.iter().enumerate() {
            if id == 2 || id == 3 || id == 6 || id == 7 {
                ux_loaded.push(out.u[i][0]);
            }
        }
        let mean: f64 = ux_loaded.iter().sum::<f64>() / ux_loaded.len() as f64;
        assert!(
            (mean - 0.01).abs() < 5e-5,
            "NLGEOM C3D8 ux={mean}, expected 0.01"
        );
        let mut sxx = 0.0;
        for s in &out.stress {
            sxx += s[0];
        }
        sxx /= out.stress.len() as f64;
        assert!(
            (sxx - 210.0).abs() < 0.5,
            "NLGEOM C3D8 sxx={sxx}, expected 210"
        );
        assert!(out.iters >= 1);
        assert!(out.residual < 1e-4, "residual={}", out.residual);
    }

    #[test]
    fn nlgeom_c3d6_small_strain_patch() {
        let inp = r#"
*HEADING
C3D6 NLGEOM patch uz=0.01
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 0, 10, 0
4, 0, 0, 10
5, 10, 0, 10
6, 0, 10, 10
*ELEMENT, TYPE=C3D6, ELSET=S
1, 1, 2, 3, 4, 5, 6
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*BOUNDARY
1, 3, 3
2, 3, 3
3, 3, 3
1, 1, 2
2, 2, 2
*STEP, NLGEOM
*STATIC
*CLOAD
4, 3, 3500
5, 3, 3500
6, 3, 3500
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let uz = out.u[out.model.node_index(4).unwrap()][2];
        assert!((uz - 0.01).abs() < 2e-5, "C3D6 NLGEOM uz={uz}");
    }

    #[test]
    fn nlgeom_c3d8_svk_stretch() {
        // nu=0, prescribed λ=1.2. E_xx=0.5*(λ²-1)=0.22, S_xx=E E_xx=22,
        // Cauchy σ_xx = λ S_xx = 26.4
        let inp = r#"
*HEADING
C3D8 SVK stretch λ=1.2
*NODE
1, 0, 0, 0
2, 1, 0, 0
3, 1, 1, 0
4, 0, 1, 0
5, 0, 0, 1
6, 1, 0, 1
7, 1, 1, 1
8, 0, 1, 1
*ELEMENT, TYPE=C3D8, ELSET=S
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
100, 0.0
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*BOUNDARY
1, 1, 3
4, 1, 1
4, 3, 3
5, 1, 2
8, 1, 1
2, 2, 3
3, 3, 3
6, 2, 2
2, 1, 1, 0.2
3, 1, 1, 0.2
6, 1, 1, 0.2
7, 1, 1, 0.2
*STEP, NLGEOM
*STATIC
*END STEP
"#;
        let out = solve_native(inp).unwrap();
        let u2 = out.u[out.model.node_index(2).unwrap()][0];
        assert!((u2 - 0.2).abs() < 1e-12, "prescribed ux={u2}");
        let sxx = out.stress[out.model.node_index(7).unwrap()][0];
        assert!(
            (sxx - 26.4).abs() < 0.05,
            "SVK Cauchy sxx={sxx}, expected 26.4"
        );
        let exx = out.strain[out.model.node_index(7).unwrap()][0];
        assert!((exx - 0.22).abs() < 1e-6, "GL Exx={exx}, expected 0.22");
        assert!(out.residual < 1e-8, "residual={}", out.residual);
    }
}
