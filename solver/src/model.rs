use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ElemKind {
    Hex8,
    Tet4,
    Quad4Ps,
    Quad4Pe,
    Tri3Ps,
    Tri3Pe,
}

impl ElemKind {
    pub fn from_ccx(name: &str) -> Option<Self> {
        Some(match name {
            "C3D8" | "C3D8R" | "C3D8I" => Self::Hex8,
            "C3D4" => Self::Tet4,
            "CPS4" | "CPS4R" | "S4" | "S4R" => Self::Quad4Ps,
            "CPE4" | "CPE4R" | "CPE4I" => Self::Quad4Pe,
            "CPS3" | "S3" => Self::Tri3Ps,
            "CPE3" => Self::Tri3Pe,
            _ => return None,
        })
    }

    pub fn nnodes(self) -> usize {
        match self {
            Self::Hex8 => 8,
            Self::Tet4 => 4,
            Self::Quad4Ps | Self::Quad4Pe => 4,
            Self::Tri3Ps | Self::Tri3Pe => 3,
        }
    }

    pub fn spatial_dim(self) -> usize {
        match self {
            Self::Hex8 | Self::Tet4 => 3,
            _ => 2,
        }
    }

    pub fn is_plane_strain(self) -> bool {
        matches!(self, Self::Quad4Pe | Self::Tri3Pe)
    }

    pub fn ccx_name(self) -> &'static str {
        match self {
            Self::Hex8 => "C3D8",
            Self::Tet4 => "C3D4",
            Self::Quad4Ps => "CPS4",
            Self::Quad4Pe => "CPE4",
            Self::Tri3Ps => "CPS3",
            Self::Tri3Pe => "CPE3",
        }
    }

    /// cgx / FRD element type numbers
    pub fn frd_type(self) -> i32 {
        match self {
            Self::Hex8 => 1,
            Self::Tet4 => 3,
            Self::Tri3Ps | Self::Tri3Pe => 7,
            Self::Quad4Ps | Self::Quad4Pe => 9,
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
}

impl Default for Material {
    fn default() -> Self {
        Self {
            e: 210000.0,
            nu: 0.3,
            density: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Boundary {
    pub node: i32,
    pub dof: usize, // 0,1,2
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
    pub nsets: HashMap<String, Vec<i32>>,
    pub elsets: HashMap<String, Vec<i32>>,
    pub bcs: Vec<Boundary>,
    pub cloads: Vec<Cload>,
    pub dloads: Vec<Dload>,
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
            nsets: HashMap::new(),
            elsets: HashMap::new(),
            bcs: Vec::new(),
            cloads: Vec::new(),
            dloads: Vec::new(),
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
}
