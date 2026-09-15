use std::collections::{HashMap, HashSet};

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
    Wedge15,
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
    Cax4,
    Cax4R,
    Cax8,
    Cax8R,
    Cax3,
    Cax6,
    Beam31,
    Beam32,
    Shell4,
    Shell4R,
    Shell3,
    Shell8,
    Shell8R,
    Shell6,
    Mem3,
    Mem4,
    Mem4R,
    Mem6,
    Mem8,
    Truss2,
    Truss3,
    SpringA,
    Tet10T,
    Mass,
    RotaryI,
    DashpotA,
    GapUni,
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
            "C3D15" => Self::Wedge15,
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
            "CAX4" => Self::Cax4,
            "CAX4R" => Self::Cax4R,
            "CAX8" => Self::Cax8,
            "CAX8R" => Self::Cax8R,
            "CAX3" => Self::Cax3,
            "CAX6" => Self::Cax6,
            "C3D10T" => Self::Tet10T,
            "C3D20RI" => Self::Hex20R,
            "MASS" => Self::Mass,
            "ROTARYI" => Self::RotaryI,
            "DASHPOTA" | "DASHPOT" => Self::DashpotA,
            "GAPUNI" => Self::GapUni,
            "S4" => Self::Shell4,
            "S4R" => Self::Shell4R,
            "S3" | "S3R" | "STRI3" => Self::Shell3,
            "S8" => Self::Shell8,
            "S8R" => Self::Shell8R,
            "S6" | "STRI65" => Self::Shell6,
            "M3D3" => Self::Mem3,
            "M3D4" => Self::Mem4,
            "M3D4R" => Self::Mem4R,
            "M3D6" => Self::Mem6,
            "M3D8" | "M3D8R" => Self::Mem8,
            "B31" | "B31R" | "B21" | "B21R" => Self::Beam31,
            "B32" | "B32R" | "B22" | "B22R" => Self::Beam32,
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
            Self::Tet10 | Self::Tet10T => 10,
            Self::Wedge6 => 6,
            Self::Wedge15 => 15,
            Self::Quad4Ps | Self::Quad4Pe | Self::Cax4 | Self::Cax4R | Self::Mem4 | Self::Mem4R => 4,
            Self::Quad8Ps
            | Self::Quad8Pe
            | Self::Quad8RPs
            | Self::Quad8RPe
            | Self::Cax8
            | Self::Cax8R
            | Self::Mem8 => 8,
            Self::Tri3Ps | Self::Tri3Pe | Self::Mem3 | Self::Cax3 => 3,
            Self::Tri6Ps | Self::Tri6Pe | Self::Shell6 | Self::Mem6 | Self::Cax6 => 6,
            Self::Beam31 => 2,
            Self::Beam32 => 3,
            Self::Shell4 | Self::Shell4R => 4,
            Self::Shell3 => 3,
            Self::Shell8 | Self::Shell8R => 8,
            Self::Truss2 | Self::SpringA | Self::DashpotA | Self::GapUni => 2,
            Self::Truss3 => 3,
            Self::Mass | Self::RotaryI => 1,
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
            | Self::Tet10T
            | Self::Wedge6
            | Self::Wedge15
            | Self::Beam31
            | Self::Beam32
            | Self::Truss2
            | Self::Truss3
            | Self::SpringA
            | Self::Mass
            | Self::RotaryI
            | Self::DashpotA
            | Self::GapUni
            | Self::Shell4
            | Self::Shell4R
            | Self::Shell3
            | Self::Shell8
            | Self::Shell8R
            | Self::Shell6
            | Self::Mem3
            | Self::Mem4
            | Self::Mem4R
            | Self::Mem6
            | Self::Mem8 => 3,
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
            Self::RotaryI => 6,
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

    pub fn is_membrane(self) -> bool {
        matches!(
            self,
            Self::Mem3 | Self::Mem4 | Self::Mem4R | Self::Mem6 | Self::Mem8
        )
    }

    pub fn is_axisym(self) -> bool {
        matches!(
            self,
            Self::Cax4 | Self::Cax4R | Self::Cax8 | Self::Cax8R | Self::Cax3 | Self::Cax6
        )
    }

    pub fn is_continuum3d(self) -> bool {
        matches!(
            self,
            Self::Hex8
                | Self::Hex8I
                | Self::Hex8R
                | Self::Hex20
                | Self::Hex20R
                | Self::Tet4
                | Self::Tet10
                | Self::Tet10T
                | Self::Wedge6
                | Self::Wedge15
        )
    }

    pub fn is_truss(self) -> bool {
        matches!(self, Self::Truss2 | Self::Truss3)
    }

    pub fn is_spring(self) -> bool {
        matches!(self, Self::SpringA)
    }

    pub fn is_point(self) -> bool {
        matches!(self, Self::Mass | Self::RotaryI)
    }

    pub fn is_dashpot(self) -> bool {
        matches!(self, Self::DashpotA)
    }

    pub fn is_gap(self) -> bool {
        matches!(self, Self::GapUni)
    }

    /// MASS / ROTARYI / DASHPOT / GAP have no continuum *MATERIAL.
    pub fn is_special(self) -> bool {
        self.is_point() || self.is_dashpot() || self.is_gap()
    }

    /// Constitutive *ELASTIC is required (not springs or concentrated mass).
    pub fn needs_material(self) -> bool {
        !self.is_special() && !self.is_spring()
    }

    pub fn is_quadratic(self) -> bool {
        matches!(
            self,
            Self::Hex20
                | Self::Hex20R
                | Self::Tet10
                | Self::Tet10T
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
                | Self::Wedge15
                | Self::Cax8
                | Self::Cax8R
                | Self::Cax6
                | Self::Mem6
                | Self::Mem8
        )
    }

    pub fn reduced_int(self) -> bool {
        matches!(
            self,
            Self::Hex20R
                | Self::Hex8R
                | Self::Quad8RPs
                | Self::Quad8RPe
                | Self::Shell4R
                | Self::Shell8R
                | Self::Cax4R
                | Self::Cax8R
                | Self::Mem4R
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
            Self::Tet10T => "C3D10T",
            Self::Wedge6 => "C3D6",
            Self::Wedge15 => "C3D15",
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
            Self::Cax4 => "CAX4",
            Self::Cax4R => "CAX4R",
            Self::Cax8 => "CAX8",
            Self::Cax8R => "CAX8R",
            Self::Cax3 => "CAX3",
            Self::Cax6 => "CAX6",
            Self::Beam31 => "B31",
            Self::Beam32 => "B32",
            Self::Shell4 => "S4",
            Self::Shell4R => "S4R",
            Self::Shell3 => "S3",
            Self::Shell8 => "S8",
            Self::Shell8R => "S8R",
            Self::Shell6 => "S6",
            Self::Mem3 => "M3D3",
            Self::Mem4 => "M3D4",
            Self::Mem4R => "M3D4R",
            Self::Mem6 => "M3D6",
            Self::Mem8 => "M3D8",
            Self::Truss2 => "T3D2",
            Self::Truss3 => "T3D3",
            Self::SpringA => "SPRINGA",
            Self::Mass => "MASS",
            Self::RotaryI => "ROTARYI",
            Self::DashpotA => "DASHPOTA",
            Self::GapUni => "GAPUNI",
        }
    }

    /// cgx / FRD element type numbers
    pub fn frd_type(self) -> i32 {
        match self {
            Self::Hex8 | Self::Hex8I | Self::Hex8R => 1,
            Self::Wedge6 => 2,
            Self::Tet4 => 3,
            Self::Hex20 | Self::Hex20R => 4,
            Self::Wedge15 => 5,
            Self::Tet10 | Self::Tet10T => 6,
            Self::Tri3Ps | Self::Tri3Pe | Self::Shell3 | Self::Mem3 | Self::Cax3 => 7,
            Self::Tri6Ps | Self::Tri6Pe | Self::Shell6 | Self::Mem6 | Self::Cax6 => 8,
            Self::Quad4Ps
            | Self::Quad4Pe
            | Self::Shell4
            | Self::Shell4R
            | Self::Cax4
            | Self::Cax4R
            | Self::Mem4
            | Self::Mem4R => 9,
            Self::Quad8Ps
            | Self::Quad8Pe
            | Self::Quad8RPs
            | Self::Quad8RPe
            | Self::Shell8
            | Self::Shell8R
            | Self::Cax8
            | Self::Cax8R
            | Self::Mem8 => 10,
            Self::Beam31
            | Self::Truss2
            | Self::SpringA
            | Self::Mass
            | Self::RotaryI
            | Self::DashpotA
            | Self::GapUni => 11,
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
    pub conductivity: f64,
    pub specific_heat: f64,
}

impl Default for Material {
    fn default() -> Self {
        Self {
            e: 210000.0,
            nu: 0.3,
            density: 0.0,
            alpha: 0.0,
            tref: 0.0,
            conductivity: 0.0,
            specific_heat: 0.0,
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

#[derive(Clone, Copy, Debug)]
pub struct GapSection {
    pub clearance: f64,
    pub k: f64,
    pub dir: [f64; 3],
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

    /// Thin-walled box: outer a×b, wall thicknesses t_bottom, t_top, t_left, t_right.
    pub fn box_sec(
        width: f64,
        height: f64,
        t_bot: f64,
        t_top: f64,
        t_left: f64,
        t_right: f64,
        n1: [f64; 3],
    ) -> Self {
        let a = width.abs().max(1e-16);
        let b = height.abs().max(1e-16);
        let tb = t_bot.abs().clamp(1e-16, b * 0.49);
        let tt = t_top.abs().clamp(1e-16, b * 0.49);
        let tl = t_left.abs().clamp(1e-16, a * 0.49);
        let tr = t_right.abs().clamp(1e-16, a * 0.49);
        let ai = (a - tl - tr).max(0.0);
        let bi = (b - tb - tt).max(0.0);
        let area = (a * b - ai * bi).max(1e-16);
        let i11 = (a * b * b * b - ai * bi * bi * bi) / 12.0;
        let i22 = (b * a * a * a - bi * ai * ai * ai) / 12.0;
        Self {
            a,
            b,
            n1,
            area,
            i11: i11.abs().max(1e-30),
            i12: 0.0,
            i22: i22.abs().max(1e-30),
            jtor: torsion_rect(a, b) - torsion_rect(ai.max(1e-16), bi.max(1e-16)),
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

#[derive(Clone, Debug)]
pub struct Cload {
    pub node: i32,
    pub dof: usize,
    pub mag: f64,
    pub amplitude: String,
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
    /// `*DLOAD, CENTRIF`: ω² and two points on the rotation axis.
    /// `elems` empty → all elements.
    Centrif {
        omega2: f64,
        p1: [f64; 3],
        p2: [f64; 3],
        elems: Vec<i32>,
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
    pub rot_node: Option<i32>,
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
    pub cylindrical: bool,
    pub origin: [f64; 3],
    pub axis: [f64; 3],
}

#[derive(Clone, Debug)]
pub struct SurfaceInteraction {
    pub kn: f64,
    pub mu: f64,
}

#[derive(Clone, Debug)]
pub struct ContactPair {
    pub slave: String,
    pub master: String,
    pub kn: f64,
    pub mu: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct ThermalBc {
    pub node: i32,
    pub value: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Cflux {
    pub node: i32,
    pub mag: f64,
}

#[derive(Clone, Copy, Debug)]
pub enum FluxKind {
    Body,
    Face(i32),
}

#[derive(Clone, Copy, Debug)]
pub struct Dflux {
    pub elem: i32,
    pub kind: FluxKind,
    pub mag: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Film {
    pub elem: i32,
    pub face: i32,
    pub t_inf: f64,
    pub h: f64,
}

#[derive(Clone, Debug)]
pub struct Amplitude {
    pub name: String,
    pub points: Vec<(f64, f64)>,
}

impl Amplitude {
    pub fn value_at(&self, t: f64) -> f64 {
        if self.points.is_empty() {
            return 1.0;
        }
        if t <= self.points[0].0 {
            return self.points[0].1;
        }
        for w in self.points.windows(2) {
            let (t0, v0) = w[0];
            let (t1, v1) = w[1];
            if t <= t1 {
                let d = t1 - t0;
                if d.abs() < 1e-18 {
                    return v1;
                }
                let a = (t - t0) / d;
                return v0 + a * (v1 - v0);
            }
        }
        self.points[self.points.len() - 1].1
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InitKind {
    Displacement,
    Velocity,
    Temperature,
}

#[derive(Clone, Copy, Debug)]
pub struct InitCond {
    pub node: i32,
    pub dof: usize,
    pub value: f64,
    pub kind: InitKind,
}

/// CalculiX `*STATIC, RIKS` data line: initial Δλ, period, min, max, max increments.
#[derive(Clone, Copy, Debug)]
pub struct RiksCtrl {
    pub dlam: f64,
    pub period: f64,
    pub dlam_min: f64,
    pub dlam_max: f64,
    pub max_inc: usize,
}

impl Default for RiksCtrl {
    fn default() -> Self {
        Self {
            dlam: 1.0,
            period: 1.0,
            dlam_min: 1e-5,
            dlam_max: 1.0,
            max_inc: 100,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Procedure {
    Static {
        nlgeom: bool,
        increments: usize,
        riks: bool,
    },
    Frequency { nmodes: usize },
    Buckle { nmodes: usize },
    HeatTransfer { steady: bool, dt: f64, period: f64 },
    Dynamic { dt: f64, period: f64 },
}

#[derive(Clone, Debug)]
pub struct AnalysisStep {
    pub procedure: Procedure,
    pub n_cload: usize,
    pub n_dload: usize,
    pub n_bc: usize,
    pub cload_from: usize,
    pub dload_from: usize,
    pub bc_from: usize,
}

impl Default for Procedure {
    fn default() -> Self {
        Self::Static {
            nlgeom: false,
            increments: 1,
            riks: false,
        }
    }
}

impl Procedure {
    pub fn name(&self) -> &'static str {
        match self {
            Self::Static {
                nlgeom: true,
                riks: true,
                ..
            } => "STATIC, NLGEOM, RIKS",
            Self::Static { riks: true, .. } => "STATIC, RIKS",
            Self::Static { nlgeom: true, .. } => "STATIC, NLGEOM",
            Self::Static { .. } => "STATIC",
            Self::Frequency { .. } => "FREQUENCY",
            Self::Buckle { .. } => "BUCKLE",
            Self::HeatTransfer { steady: true, .. } => "HEAT TRANSFER, STEADY STATE",
            Self::HeatTransfer { .. } => "HEAT TRANSFER",
            Self::Dynamic { .. } => "DYNAMIC",
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
    pub elset_mass: HashMap<String, f64>,
    pub elset_rotary: HashMap<String, [f64; 6]>,
    pub elset_dashpot: HashMap<String, f64>,
    pub elset_gap: HashMap<String, GapSection>,
    pub temperatures: HashMap<i32, f64>,
    pub plastic: HashMap<String, Vec<(f64, f64)>>, // material -> [(peeq, sy)]
    pub interactions: HashMap<String, SurfaceInteraction>,
    pub contact_pairs: Vec<ContactPair>,
    pub thermal_bcs: Vec<ThermalBc>,
    pub cfluxes: Vec<Cflux>,
    pub dfluxes: Vec<Dflux>,
    pub films: Vec<Film>,
    pub amplitudes: Vec<Amplitude>,
    pub init: Vec<InitCond>,
    pub damp_alpha: f64,
    pub damp_beta: f64,
    pub procedure: Procedure,
    pub steps: Vec<AnalysisStep>,
    pub max_newton: usize,
    pub newton_tol: f64,
    pub riks: Option<RiksCtrl>,
    pub dim: usize,
    pub output_u: bool,
    pub output_s: bool,
    pub output_rf: bool,
    pub output_e: bool,
    pub output_nt: bool,
    pub warnings: Vec<String>,
    /// Displacements at the start of the current step (multi-step NLGEOM).
    pub u_start: Vec<[f64; 3]>,
    /// Nodal load vector at the start of the current step (ndof = 3*nnode).
    pub f_start: Vec<f64>,
    /// `*STEP, INC=` cap on static increments (CalculiX default 100).
    pub max_inc: usize,
    pub static_dt: f64,
    pub static_period: f64,
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
            elset_mass: HashMap::new(),
            elset_rotary: HashMap::new(),
            elset_dashpot: HashMap::new(),
            elset_gap: HashMap::new(),
            temperatures: HashMap::new(),
            plastic: HashMap::new(),
            interactions: HashMap::new(),
            contact_pairs: Vec::new(),
            thermal_bcs: Vec::new(),
            cfluxes: Vec::new(),
            dfluxes: Vec::new(),
            films: Vec::new(),
            amplitudes: Vec::new(),
            init: Vec::new(),
            damp_alpha: 0.0,
            damp_beta: 0.0,
            procedure: Procedure::default(),
            steps: Vec::new(),
            max_newton: 25,
            newton_tol: 1e-8,
            riks: None,
            dim: 3,
            output_u: true,
            output_s: true,
            output_rf: true,
            output_e: false,
            output_nt: false,
            warnings: Vec::new(),
            u_start: Vec::new(),
            f_start: Vec::new(),
            max_inc: 100,
            static_dt: 1.0,
            static_period: 1.0,
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
            || self.elements.iter().any(|e| e.kind == ElemKind::RotaryI)
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
        self.apply_named_elset_sections();
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

    /// CalculiX/Mecway: `*ELEMENT` often has no `ELSET=`; properties come from a
    /// later `*ELSET` plus `*SHELL SECTION` / `*MASS` / `*BEAM SECTION` on that set.
    /// Rebind each element to a sectioned elset it belongs to.
    fn apply_named_elset_sections(&mut self) {
        let mut sectioned: HashSet<String> = HashSet::new();
        sectioned.extend(self.elset_material.keys().cloned());
        sectioned.extend(self.elset_thickness.keys().cloned());
        sectioned.extend(self.elset_beam.keys().cloned());
        sectioned.extend(self.elset_spring.keys().cloned());
        sectioned.extend(self.elset_mass.keys().cloned());
        sectioned.extend(self.elset_rotary.keys().cloned());
        sectioned.extend(self.elset_dashpot.keys().cloned());
        sectioned.extend(self.elset_gap.keys().cloned());
        if sectioned.is_empty() {
            return;
        }
        let mut names: Vec<String> = sectioned.iter().cloned().collect();
        names.sort_by(|a, b| {
            let am = self.elset_material.contains_key(a) as u8;
            let bm = self.elset_material.contains_key(b) as u8;
            bm.cmp(&am).then(a.cmp(b))
        });
        let mut id_to_set: HashMap<i32, String> = HashMap::new();
        for name in &names {
            if let Some(ids) = self.elsets.get(name) {
                for &id in ids {
                    id_to_set.entry(id).or_insert_with(|| name.clone());
                }
            }
        }
        for el in &mut self.elements {
            if sectioned.contains(&el.elset) {
                continue;
            }
            if let Some(name) = id_to_set.get(&el.id) {
                el.elset = name.clone();
            }
        }
    }

    fn expand_transforms(&mut self) {
        let mut map = HashMap::new();
        let nsets = self.nsets.clone();
        let id_to_idx: HashMap<i32, usize> = self
            .node_ids
            .iter()
            .enumerate()
            .map(|(i, id)| (*id, i))
            .collect();
        for t in &self.transforms {
            if let Some(nodes) = nsets.get(&t.nset) {
                for &id in nodes {
                    if t.cylindrical {
                        if let Some(&i) = id_to_idx.get(&id) {
                            let p = self.coords[i];
                            map.insert(id, cylindrical_axes(t.origin, t.axis, p));
                            continue;
                        }
                    }
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
        if let Some(v) = self.nsets.get(name) {
            return Ok(v.clone());
        }
        let up = name.to_ascii_uppercase();
        if let Some(v) = self.nsets.get(&up) {
            return Ok(v.clone());
        }
        if up == "NALL" || up == "ALL" {
            return Ok(self.node_ids.clone());
        }
        crate::error::err(format!("Unbekanntes Knotenset {name}"))
    }

    pub fn expand_elset(&self, name: &str) -> crate::error::Result<Vec<i32>> {
        if let Ok(id) = name.parse::<i32>() {
            return Ok(vec![id]);
        }
        if let Some(v) = self.elsets.get(name) {
            return Ok(v.clone());
        }
        let up = name.to_ascii_uppercase();
        if let Some(v) = self.elsets.get(&up) {
            return Ok(v.clone());
        }
        if up == "EALL" || up == "ALL" {
            return Ok(self.elements.iter().map(|e| e.id).collect());
        }
        crate::error::err(format!("Unbekanntes Elementset {name}"))
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

    pub fn mass_for(&self, el: &Element) -> crate::error::Result<f64> {
        if let Some(&m) = self.elset_mass.get(&el.elset) {
            return Ok(m);
        }
        if self.elset_mass.len() == 1 {
            return Ok(*self.elset_mass.values().next().unwrap());
        }
        crate::error::err(format!(
            "Keine *MASS für Element {} (ELSET={})",
            el.id, el.elset
        ))
    }

    pub fn rotary_for(&self, el: &Element) -> crate::error::Result<[f64; 6]> {
        if let Some(&i) = self.elset_rotary.get(&el.elset) {
            return Ok(i);
        }
        if self.elset_rotary.len() == 1 {
            return Ok(*self.elset_rotary.values().next().unwrap());
        }
        crate::error::err(format!(
            "Keine *ROTARY INERTIA für Element {} (ELSET={})",
            el.id, el.elset
        ))
    }

    pub fn dashpot_for(&self, el: &Element) -> crate::error::Result<f64> {
        if let Some(&c) = self.elset_dashpot.get(&el.elset) {
            return Ok(c);
        }
        if self.elset_dashpot.len() == 1 {
            return Ok(*self.elset_dashpot.values().next().unwrap());
        }
        crate::error::err(format!(
            "Keine *DASHPOT-Dämpfung für Element {} (ELSET={})",
            el.id, el.elset
        ))
    }

    pub fn gap_for(&self, el: &Element) -> crate::error::Result<GapSection> {
        if let Some(&g) = self.elset_gap.get(&el.elset) {
            return Ok(g);
        }
        if self.elset_gap.len() == 1 {
            return Ok(*self.elset_gap.values().next().unwrap());
        }
        Ok(GapSection {
            clearance: 0.0,
            k: 1.0e8,
            dir: [0.0, 0.0, 0.0],
        })
    }

    pub fn temperature_at(&self, node: i32) -> f64 {
        self.temperatures.get(&node).copied().unwrap_or(0.0)
    }

    pub fn plastic_for(&self, el: &Element) -> Option<&[(f64, f64)]> {
        if let Some(mname) = self.elset_material.get(&el.elset) {
            if let Some(p) = self.plastic.get(mname) {
                return Some(p.as_slice());
            }
        }
        if self.plastic.len() == 1 {
            return self.plastic.values().next().map(|v| v.as_slice());
        }
        None
    }

    pub fn has_plastic(&self) -> bool {
        !self.plastic.is_empty()
    }

    pub fn has_contact(&self) -> bool {
        !self.contact_pairs.is_empty()
    }

    pub fn amp_value(&self, name: &str, t: f64) -> f64 {
        if name.is_empty() {
            return 1.0;
        }
        let n = name.to_ascii_uppercase();
        self.amplitudes
            .iter()
            .find(|a| a.name == n)
            .map(|a| a.value_at(t))
            .unwrap_or(1.0)
    }

    pub fn find_amplitude(&self, name: &str) -> usize {
        let n = name.to_ascii_uppercase();
        self.amplitudes
            .iter()
            .position(|a| a.name == n)
            .map(|i| i + 1)
            .unwrap_or(0)
    }
}

fn cylindrical_axes(origin: [f64; 3], axis: [f64; 3], p: [f64; 3]) -> [[f64; 3]; 3] {
    let al = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let e3 = if al < 1e-18 {
        [0.0, 0.0, 1.0]
    } else {
        [axis[0] / al, axis[1] / al, axis[2] / al]
    };
    let r = [p[0] - origin[0], p[1] - origin[1], p[2] - origin[2]];
    let proj = r[0] * e3[0] + r[1] * e3[1] + r[2] * e3[2];
    let mut e1 = [r[0] - proj * e3[0], r[1] - proj * e3[1], r[2] - proj * e3[2]];
    let n1 = (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]).sqrt();
    if n1 < 1e-12 {
        let helper = if e3[2].abs() < 0.9 {
            [0.0, 0.0, 1.0]
        } else {
            [1.0, 0.0, 0.0]
        };
        e1 = [
            helper[1] * e3[2] - helper[2] * e3[1],
            helper[2] * e3[0] - helper[0] * e3[2],
            helper[0] * e3[1] - helper[1] * e3[0],
        ];
        let n = (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]).sqrt().max(1e-30);
        e1 = [e1[0] / n, e1[1] / n, e1[2] / n];
    } else {
        e1 = [e1[0] / n1, e1[1] / n1, e1[2] / n1];
    }
    let e2 = [
        e3[1] * e1[2] - e3[2] * e1[1],
        e3[2] * e1[0] - e3[0] * e1[2],
        e3[0] * e1[1] - e3[1] * e1[0],
    ];
    [
        [e1[0], e2[0], e3[0]],
        [e1[1], e2[1], e3[1]],
        [e1[2], e2[2], e3[2]],
    ]
}
