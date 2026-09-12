use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ElemKind {
    Hex8,
    Hex8I,
    Hex8R,
    Hex20,
    Hex20R,
    Tet4,
    Tet10,
    Wedge6,
    Quad4Ps,
    Quad4Pe,
    Quad8Ps,
    Quad8Pe,
    Quad8RPs,
    Quad8RPe,
    Tri3Ps,
    Tri3Pe,
    Tri6Ps,
    Tri6Pe,
    Beam31,
    Beam32,
    Shell4,
    Shell4R,
    Shell3,
    Shell8,
    Shell8R,
    Shell6,
    Truss2,
    Truss3,
    SpringA,
}

impl ElemKind {
    pub fn from_ccx(name: &str) -> Option<Self> {
        Some(match name {
            "C3D8" => Self::Hex8,
            "C3D8R" => Self::Hex8R,
            "C3D8I" => Self::Hex8I,
            "C3D20" => Self::Hex20,
            "C3D20R" => Self::Hex20R,
            "C3D4" => Self::Tet4,
            "C3D10" | "C3D10M" => Self::Tet10,
            "C3D6" => Self::Wedge6,
            "CPS4" | "CPS4R" => Self::Quad4Ps,
            "CPE4" | "CPE4R" | "CPE4I" => Self::Quad4Pe,
            "CPS8" => Self::Quad8Ps,
            "CPE8" => Self::Quad8Pe,
            "CPS8R" => Self::Quad8RPs,
            "CPE8R" => Self::Quad8RPe,
            "CPS3" => Self::Tri3Ps,
            "CPE3" => Self::Tri3Pe,
            "CPS6" => Self::Tri6Ps,
            "CPE6" => Self::Tri6Pe,
            "S4" => Self::Shell4,
            "S4R" => Self::Shell4R,
            "S3" | "S3R" | "STRI3" => Self::Shell3,
            "S8" => Self::Shell8,
            "S8R" => Self::Shell8R,
            "S6" | "STRI65" => Self::Shell6,
            "B31" | "B31R" => Self::Beam31,
            "B32" | "B32R" => Self::Beam32,
            "T3D2" | "T2D2" => Self::Truss2,
            "T3D3" => Self::Truss3,
            "SPRINGA" | "SPRING2" => Self::SpringA,
            _ => return None,
        })
    }

    pub fn nnodes(self) -> usize {
        match self {
            Self::Hex8 | Self::Hex8I | Self::Hex8R => 8,
            Self::Hex20 | Self::Hex20R => 20,
            Self::Tet4 => 4,
            Self::Tet10 => 10,
            Self::Wedge6 => 6,
            Self::Quad4Ps | Self::Quad4Pe => 4,
            Self::Quad8Ps | Self::Quad8Pe | Self::Quad8RPs | Self::Quad8RPe => 8,
            Self::Tri3Ps | Self::Tri3Pe => 3,
            Self::Tri6Ps | Self::Tri6Pe | Self::Shell6 => 6,
            Self::Beam31 => 2,
            Self::Beam32 => 3,
            Self::Shell4 | Self::Shell4R => 4,
            Self::Shell3 => 3,
            Self::Shell8 | Self::Shell8R => 8,
            Self::Truss2 | Self::SpringA => 2,
            Self::Truss3 => 3,
        }
    }

    pub fn spatial_dim(self) -> usize {
        match self {
            Self::Hex8
            | Self::Hex8I
            | Self::Hex8R
            | Self::Hex20
            | Self::Hex20R
            | Self::Tet4
            | Self::Tet10
            | Self::Wedge6
            | Self::Beam31
            | Self::Beam32
            | Self::Truss2
            | Self::Truss3
            | Self::SpringA
            | Self::Shell4
            | Self::Shell4R
            | Self::Shell3
            | Self::Shell8
            | Self::Shell8R
            | Self::Shell6 => 3,
            _ => 2,
        }
    }

    pub fn ndof_per_node(self) -> usize {
        match self {
            Self::Beam31
            | Self::Beam32
            | Self::Shell4
            | Self::Shell4R
            | Self::Shell3
            | Self::Shell8
            | Self::Shell8R
            | Self::Shell6 => 6,
            k if k.spatial_dim() == 3 => 3,
            _ => 2,
        }
    }

    pub fn is_beam(self) -> bool {
        matches!(self, Self::Beam31 | Self::Beam32)
    }

    pub fn is_shell(self) -> bool {
        matches!(
            self,
            Self::Shell4
                | Self::Shell4R
                | Self::Shell3
                | Self::Shell8
                | Self::Shell8R
                | Self::Shell6
        )
    }

    pub fn is_truss(self) -> bool {
        matches!(self, Self::Truss2 | Self::Truss3)
    }

    pub fn is_spring(self) -> bool {
        matches!(self, Self::SpringA)
    }

    pub fn is_quadratic(self) -> bool {
        matches!(
            self,
            Self::Hex20
                | Self::Hex20R
                | Self::Tet10
                | Self::Quad8Ps
                | Self::Quad8Pe
                | Self::Quad8RPs
                | Self::Quad8RPe
                | Self::Tri6Ps
                | Self::Tri6Pe
                | Self::Beam32
                | Self::Shell8
                | Self::Shell8R
                | Self::Shell6
        )
    }

    pub fn reduced_int(self) -> bool {
        matches!(
            self,
            Self::Hex20R | Self::Hex8R | Self::Quad8RPs | Self::Quad8RPe | Self::Shell4R | Self::Shell8R
        )
    }

    pub fn is_plane_strain(self) -> bool {
        matches!(
            self,
            Self::Quad4Pe | Self::Quad8Pe | Self::Quad8RPe | Self::Tri3Pe | Self::Tri6Pe
        )
    }

    pub fn ccx_name(self) -> &'static str {
        match self {
            Self::Hex8 => "C3D8",
            Self::Hex8I => "C3D8I",
            Self::Hex8R => "C3D8R",
            Self::Hex20 => "C3D20",
            Self::Hex20R => "C3D20R",
            Self::Tet4 => "C3D4",
            Self::Tet10 => "C3D10",
            Self::Wedge6 => "C3D6",
            Self::Quad4Ps => "CPS4",
            Self::Quad4Pe => "CPE4",
            Self::Quad8Ps => "CPS8",
            Self::Quad8Pe => "CPE8",
            Self::Quad8RPs => "CPS8R",
            Self::Quad8RPe => "CPE8R",
            Self::Tri3Ps => "CPS3",
            Self::Tri3Pe => "CPE3",
            Self::Tri6Ps => "CPS6",
            Self::Tri6Pe => "CPE6",
            Self::Beam31 => "B31",
            Self::Beam32 => "B32",
            Self::Shell4 => "S4",
            Self::Shell4R => "S4R",
            Self::Shell3 => "S3",
            Self::Shell8 => "S8",
            Self::Shell8R => "S8R",
            Self::Shell6 => "S6",
            Self::Truss2 => "T3D2",
            Self::Truss3 => "T3D3",
            Self::SpringA => "SPRINGA",
        }
    }

    /// cgx / FRD element type numbers
    pub fn frd_type(self) -> i32 {
        match self {
            Self::Hex8 | Self::Hex8I | Self::Hex8R => 1,
            Self::Wedge6 => 2,
            Self::Hex20 | Self::Hex20R => 4,
            Self::Tet4 => 3,
            Self::Tet10 => 6,
            Self::Tri3Ps | Self::Tri3Pe | Self::Shell3 => 7,
            Self::Tri6Ps | Self::Tri6Pe | Self::Shell6 => 8,
            Self::Quad4Ps | Self::Quad4Pe | Self::Shell4 | Self::Shell4R => 9,
            Self::Quad8Ps
            | Self::Quad8Pe
            | Self::Quad8RPs
            | Self::Quad8RPe
            | Self::Shell8
            | Self::Shell8R => 10,
            Self::Beam31 | Self::Truss2 | Self::SpringA => 11,
            Self::Beam32 | Self::Truss3 => 12,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Element {
    pub id: i32,
    pub kind: ElemKind,
    pub nodes: Vec<i32>,
    pub elset: String,
}

#[derive(Clone, Copy, Debug)]
pub struct Material {
    pub e: f64,
    pub nu: f64,
    pub density: f64,
    pub alpha: f64,
    pub tref: f64,
}

impl Default for Material {
    fn default() -> Self {
        Self {
            e: 210000.0,
            nu: 0.3,
            density: 0.0,
            alpha: 0.0,
            tref: 0.0,
        }
    }
}

/// Abaqus/CalculiX beam section. Local 1 = tangent t, 2 = n1, 3 = n2 = t × n1.
#[derive(Clone, Copy, Debug)]
pub struct BeamSection {
    pub a: f64,
    pub b: f64,
    pub n1: [f64; 3],
    pub area: f64,
    pub i11: f64,
    pub i12: f64,
    pub i22: f64,
    pub jtor: f64,
    pub k11: f64,
    pub k22: f64,
}

impl BeamSection {
    pub fn rect(width_n1: f64, height_n2: f64, n1: [f64; 3]) -> Self {
        let a = width_n1.abs().max(1e-16);
        let b = height_n2.abs().max(1e-16);
        Self {
            a,
            b,
            n1,
            area: a * b,
            i11: a * b * b * b / 12.0,
            i12: 0.0,
            i22: b * a * a * a / 12.0,
            jtor: torsion_rect(a, b),
            k11: 5.0 / 6.0,
            k22: 5.0 / 6.0,
        }
    }

    pub fn circ(radius: f64, n1: [f64; 3]) -> Self {
        let r = radius.abs().max(1e-16);
        let r2 = r * r;
        let r4 = r2 * r2;
        Self {
            a: r,
            b: r,
            n1,
            area: std::f64::consts::PI * r2,
            i11: std::f64::consts::PI * r4 / 4.0,
            i12: 0.0,
            i22: std::f64::consts::PI * r4 / 4.0,
            jtor: std::f64::consts::PI * r4 / 2.0,
            k11: 0.9,
            k22: 0.9,
        }
    }

    pub fn pipe(r_outer: f64, thickness: f64, n1: [f64; 3]) -> Self {
        let ro = r_outer.abs().max(1e-16);
        let t = thickness.abs().clamp(1e-16, ro);
        let ri = (ro - t).max(0.0);
        let ro2 = ro * ro;
        let ri2 = ri * ri;
        let ro4 = ro2 * ro2;
        let ri4 = ri2 * ri2;
        Self {
            a: ro,
            b: ro,
            n1,
            area: std::f64::consts::PI * (ro2 - ri2),
            i11: std::f64::consts::PI / 4.0 * (ro4 - ri4),
            i12: 0.0,
            i22: std::f64::consts::PI / 4.0 * (ro4 - ri4),
            jtor: std::f64::consts::PI / 2.0 * (ro4 - ri4),
            k11: 0.5,
            k22: 0.5,
        }
    }

    pub fn general(area: f64, i11: f64, i12: f64, i22: f64, jtor: f64, n1: [f64; 3]) -> Self {
        let a = area.abs().max(1e-16);
        // Viewer fallback: square of equal area
        let side = a.sqrt();
        Self {
            a: side,
            b: side,
            n1,
            area: a,
            i11: i11.abs().max(1e-30),
            i12,
            i22: i22.abs().max(1e-30),
            jtor: jtor.abs().max(1e-30),
            k11: 5.0 / 6.0,
            k22: 5.0 / 6.0,
        }
    }
}

fn torsion_rect(a: f64, b: f64) -> f64 {
    let (aa, bb) = if a >= b { (a, b) } else { (b, a) };
    let ratio = bb / aa;
    let ratio4 = ratio * ratio * ratio * ratio;
    aa * bb * bb * bb * (1.0 / 3.0 - 0.21 * ratio * (1.0 - ratio4 / 12.0))
}

#[derive(Clone, Copy, Debug)]
pub struct Boundary {
    pub node: i32,
    pub dof: usize, // 0..5  (u1,u2,u3,ur1,ur2,ur3)
    pub value: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Cload {
    pub node: i32,
    pub dof: usize,
    pub mag: f64,
}

#[derive(Clone, Debug)]
pub enum Dload {
    /// Face/edge pressure. face is 1-based P1..P6. `None` means body/face P for 2D surface (unused).
    Pressure {
        elem: i32,
        face: i32,
        mag: f64,
    },
    Grav {
        mag: f64,
        dir: [f64; 3],
    },
    /// Force per unit length in a global direction (Abaqus PX/PY/PZ on beams).
    BeamGlobal {
        elem: i32,
        dir: [f64; 3],
        mag: f64,
    },
}

#[derive(Clone, Debug)]
pub struct Equation {
    pub terms: Vec<(i32, usize, f64)>,
    pub rhs: f64,
}

#[derive(Clone, Debug)]
pub struct Surface {
    pub name: String,
    pub nodes: Vec<i32>,
    pub faces: Vec<(i32, i32)>,
}

#[derive(Clone, Debug)]
pub struct Tie {
    pub slave: String,
    pub master: String,
    pub position_tol: f64,
}

#[derive(Clone, Debug)]
pub struct RigidBody {
    pub nset: String,
    pub ref_node: i32,
}

#[derive(Clone, Debug)]
pub struct Coupling {
    pub ref_node: i32,
    pub surface: String,
    pub kinematic: bool,
    pub dofs: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct Transform {
    pub nset: String,
    pub axes: [[f64; 3]; 3],
}

#[derive(Clone, Debug, PartialEq)]
pub enum Procedure {
    Static { nlgeom: bool, increments: usize },
    Frequency { nmodes: usize },
    Buckle { nmodes: usize },
}

impl Default for Procedure {
    fn default() -> Self {
        Self::Static {
            nlgeom: false,
            increments: 1,
        }
    }
}

impl Procedure {
    pub fn name(&self) -> &'static str {
        match self {
            Self::Static { nlgeom: true, .. } => "STATIC, NLGEOM",
            Self::Static { .. } => "STATIC",
            Self::Frequency { .. } => "FREQUENCY",
            Self::Buckle { .. } => "BUCKLE",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Model {
    pub heading: String,
    pub node_ids: Vec<i32>,
    pub coords: Vec<[f64; 3]>,
    pub id_to_index: HashMap<i32, usize>,
    pub elements: Vec<Element>,
    pub materials: HashMap<String, Material>,
    pub elset_material: HashMap<String, String>,
    pub elset_thickness: HashMap<String, f64>,
    pub elset_beam: HashMap<String, BeamSection>,
    pub nsets: HashMap<String, Vec<i32>>,
    pub elsets: HashMap<String, Vec<i32>>,
    pub bcs: Vec<Boundary>,
    pub cloads: Vec<Cload>,
    pub dloads: Vec<Dload>,
    pub equations: Vec<Equation>,
    pub surfaces: HashMap<String, Surface>,
    pub ties: Vec<Tie>,
    pub rigid_bodies: Vec<RigidBody>,
    pub couplings: Vec<Coupling>,
    pub transforms: Vec<Transform>,
    pub node_transform: HashMap<i32, [[f64; 3]; 3]>,
    pub elset_spring: HashMap<String, f64>,
    pub temperatures: HashMap<i32, f64>,
    pub procedure: Procedure,
    pub dim: usize,
    pub output_u: bool,
    pub output_s: bool,
    pub output_rf: bool,
    pub output_e: bool,
    pub warnings: Vec<String>,
}

impl Model {
    pub fn new() -> Self {
        Self {
            heading: "Axia".into(),
            node_ids: Vec::new(),
            coords: Vec::new(),
            id_to_index: HashMap::new(),
            elements: Vec::new(),
            materials: HashMap::new(),
            elset_material: HashMap::new(),
            elset_thickness: HashMap::new(),
            elset_beam: HashMap::new(),
            nsets: HashMap::new(),
            elsets: HashMap::new(),
            bcs: Vec::new(),
            cloads: Vec::new(),
            dloads: Vec::new(),
            equations: Vec::new(),
            surfaces: HashMap::new(),
            ties: Vec::new(),
            rigid_bodies: Vec::new(),
            couplings: Vec::new(),
            transforms: Vec::new(),
            node_transform: HashMap::new(),
            elset_spring: HashMap::new(),
            temperatures: HashMap::new(),
            procedure: Procedure::default(),
            dim: 3,
            output_u: true,
            output_s: true,
            output_rf: true,
            output_e: false,
            warnings: Vec::new(),
        }
    }

    pub fn warn(&mut self, msg: impl Into<String>) {
        self.warnings.push(msg.into());
    }

    pub fn has_beams(&self) -> bool {
        self.elements.iter().any(|e| e.kind.is_beam())
    }

    pub fn has_shells(&self) -> bool {
        self.elements.iter().any(|e| e.kind.is_shell())
    }

    pub fn ndof_node(&self) -> usize {
        if self.has_beams()
            || self.has_shells()
            || !self.rigid_bodies.is_empty()
            || self.couplings.iter().any(|c| c.kinematic)
        {
            6
        } else {
            self.dim
        }
    }

    pub fn compact(&mut self) {
        let mut pairs: Vec<(i32, [f64; 3])> = self
            .node_ids
            .iter()
            .copied()
            .zip(self.coords.iter().copied())
            .collect();
        pairs.sort_by_key(|(id, _)| *id);
        pairs.dedup_by_key(|(id, _)| *id);
        self.node_ids = pairs.iter().map(|p| p.0).collect();
        self.coords = pairs.iter().map(|p| p.1).collect();
        self.id_to_index = self
            .node_ids
            .iter()
            .enumerate()
            .map(|(i, id)| (*id, i))
            .collect();
        self.elements.sort_by_key(|e| e.id);
        let all_nodes = self.node_ids.clone();
        let all_elems: Vec<i32> = self.elements.iter().map(|e| e.id).collect();
        self.nsets.insert("NALL".into(), all_nodes);
        self.elsets.insert("EALL".into(), all_elems);
        let has_2d = self.elements.iter().any(|e| e.kind.spatial_dim() == 2);
        let has_3d = self.elements.iter().any(|e| e.kind.spatial_dim() == 3);
        self.dim = if has_2d && !has_3d { 2 } else { 3 };
        if has_2d && has_3d {
            self.warn("2D- und 3D-Elemente gemischt — 3D-Freiheitsgrade werden verwendet.");
            self.dim = 3;
        }
        if self.has_beams() || self.has_shells() || !self.rigid_bodies.is_empty() {
            self.dim = 3;
        }
        self.expand_transforms();
        self.expand_surfaces();
    }

    fn expand_transforms(&mut self) {
        let mut map = HashMap::new();
        let nsets = self.nsets.clone();
        for t in &self.transforms {
            if let Some(nodes) = nsets.get(&t.nset) {
                for &id in nodes {
                    map.insert(id, t.axes);
                }
            }
        }
        self.node_transform = map;
    }

    fn expand_surfaces(&mut self) {
        let nsets = self.nsets.clone();
        let elsets = self.elsets.clone();
        for s in self.surfaces.values_mut() {
            let mut extra_nodes = Vec::new();
            for n in s.nodes.clone() {
                extra_nodes.push(n);
            }
            // faces already stored; collect nodes from named sets left as-is
            let _ = (&nsets, &elsets, extra_nodes);
        }
    }

    pub fn node_index(&self, id: i32) -> crate::error::Result<usize> {
        self.id_to_index
            .get(&id)
            .copied()
            .ok_or_else(|| crate::error::FemError(format!("Unbekannter Knoten {id}")))
    }

    pub fn expand_nset(&self, name: &str) -> crate::error::Result<Vec<i32>> {
        if let Ok(id) = name.parse::<i32>() {
            return Ok(vec![id]);
        }
        self.nsets.get(name).cloned().ok_or_else(|| {
            crate::error::FemError(format!("Unbekanntes Knotenset {name}"))
        })
    }

    pub fn expand_elset(&self, name: &str) -> crate::error::Result<Vec<i32>> {
        if let Ok(id) = name.parse::<i32>() {
            return Ok(vec![id]);
        }
        self.elsets.get(name).cloned().ok_or_else(|| {
            crate::error::FemError(format!("Unbekanntes Elementset {name}"))
        })
    }

    pub fn material_for(&self, el: &Element) -> crate::error::Result<Material> {
        if let Some(mname) = self.elset_material.get(&el.elset) {
            if let Some(m) = self.materials.get(mname) {
                return Ok(*m);
            }
        }
        if self.materials.len() == 1 {
            return Ok(*self.materials.values().next().unwrap());
        }
        if let Some(m) = self.materials.get("STEEL") {
            return Ok(*m);
        }
        crate::error::err(format!(
            "Kein Material für Element {} (ELSET={})",
            el.id, el.elset
        ))
    }

    pub fn thickness_for(&self, el: &Element) -> f64 {
        self.elset_thickness
            .get(&el.elset)
            .copied()
            .or_else(|| self.elset_thickness.values().copied().next())
            .unwrap_or(1.0)
    }

    pub fn beam_section_for(&self, el: &Element) -> crate::error::Result<BeamSection> {
        if let Some(s) = self.elset_beam.get(&el.elset) {
            return Ok(*s);
        }
        if self.elset_beam.len() == 1 {
            return Ok(*self.elset_beam.values().next().unwrap());
        }
        crate::error::err(format!(
            "Keine *BEAM SECTION für Element {} (ELSET={})",
            el.id, el.elset
        ))
    }

    pub fn spring_k_for(&self, el: &Element) -> crate::error::Result<f64> {
        if let Some(&k) = self.elset_spring.get(&el.elset) {
            return Ok(k);
        }
        if self.elset_spring.len() == 1 {
            return Ok(*self.elset_spring.values().next().unwrap());
        }
        crate::error::err(format!(
            "Keine *SPRING-Steifigkeit für Element {} (ELSET={})",
            el.id, el.elset
        ))
    }

    pub fn temperature_at(&self, node: i32) -> f64 {
        self.temperatures.get(&node).copied().unwrap_or(0.0)
    }
}
