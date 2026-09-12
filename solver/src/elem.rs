use crate::beam;
use crate::error::{err, Result};
use crate::model::{BeamSection, ElemKind};
use crate::quadratic;
use crate::shell;

pub(crate) const G2: f64 = 0.5773502691896257; // 1/sqrt(3)

pub fn invert3(a: [[f64; 3]; 3]) -> Result<([[f64; 3]; 3], f64)> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-18 {
        return err("Jakobimatrix ist singulär (entartetes Element).");
    }
    let inv_det = 1.0 / det;
    let mut b = [[0.0; 3]; 3];
    b[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * inv_det;
    b[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * inv_det;
    b[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inv_det;
    b[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * inv_det;
    b[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inv_det;
    b[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * inv_det;
    b[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * inv_det;
    b[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * inv_det;
    b[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inv_det;
    Ok((b, det))
}

pub fn invert2(a: [[f64; 2]; 2]) -> Result<([[f64; 2]; 2], f64)> {
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if det.abs() < 1e-18 {
        return err("Jakobimatrix ist singulär (entartetes Element).");
    }
    let i = 1.0 / det;
    Ok(([[a[1][1] * i, -a[0][1] * i], [-a[1][0] * i, a[0][0] * i]], det))
}

/// Isotropic 3D elasticity, Voigt [xx,yy,zz,xy,yz,zx], engineering shear.
pub fn d_iso_3d(e: f64, nu: f64) -> Result<[f64; 36]> {
    if nu >= 0.5 || nu <= -1.0 {
        return err(format!("Ungültige Querkontraktion nu={nu}"));
    }
    if e <= 0.0 {
        return err(format!("Ungültiger E-Modul E={e}"));
    }
    let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let mut d = [0.0; 36];
    let diag = lam + 2.0 * mu;
    for i in 0..3 {
        for j in 0..3 {
            d[i * 6 + j] = if i == j { diag } else { lam };
        }
        d[(i + 3) * 6 + (i + 3)] = mu;
    }
    Ok(d)
}

pub fn d_plane_stress(e: f64, nu: f64) -> Result<[f64; 9]> {
    if nu >= 1.0 || nu <= -1.0 {
        return err(format!("Ungültige Querkontraktion nu={nu}"));
    }
    if e <= 0.0 {
        return err(format!("Ungültiger E-Modul E={e}"));
    }
    let c = e / (1.0 - nu * nu);
    let mut d = [0.0; 9];
    d[0] = c;
    d[1] = c * nu;
    d[3] = c * nu;
    d[4] = c;
    d[8] = c * (1.0 - nu) / 2.0;
    Ok(d)
}

pub fn d_plane_strain(e: f64, nu: f64) -> Result<[f64; 9]> {
    let d3 = d_iso_3d(e, nu)?;
    // Extract in-plane block: xx,yy,xy from 6x6
    let mut d = [0.0; 9];
    d[0] = d3[0]; // xx-xx
    d[1] = d3[1]; // xx-yy
    d[3] = d3[6]; // yy-xx
    d[4] = d3[7]; // yy-yy
    d[8] = d3[21]; // xy-xy (index 3,3 in 6x6 = 3*6+3 = 21)
    Ok(d)
}

pub(crate) fn gemm_bt_d_b(ke: &mut [f64], n: usize, b: &[f64], nrow: usize, d: &[f64], w: f64) {
    // tmp = D * B  (nrow x n)
    let mut tmp = vec![0.0; nrow * n];
    for i in 0..nrow {
        for j in 0..n {
            let mut s = 0.0;
            for k in 0..nrow {
                s += d[i * nrow + k] * b[k * n + j];
            }
            tmp[i * n + j] = s;
        }
    }
    // ke += w * B^T * tmp
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in 0..nrow {
                s += b[k * n + i] * tmp[k * n + j];
            }
            ke[i * n + j] += w * s;
        }
    }
}

pub(crate) fn fill_b3(b: &mut [f64], nnode: usize, dndx: &[[f64; 3]]) {
    // B 6 x (3*nnode), row-major
    let n = 3 * nnode;
    for i in 0..nnode {
        let c = 3 * i;
        let dx = dndx[i][0];
        let dy = dndx[i][1];
        let dz = dndx[i][2];
        b[0 * n + c] = dx;
        b[1 * n + c + 1] = dy;
        b[2 * n + c + 2] = dz;
        b[3 * n + c] = dy;
        b[3 * n + c + 1] = dx;
        b[4 * n + c + 1] = dz;
        b[4 * n + c + 2] = dy;
        b[5 * n + c] = dz;
        b[5 * n + c + 2] = dx;
    }
}

pub(crate) fn fill_b2(b: &mut [f64], nnode: usize, dndx: &[[f64; 2]]) {
    let n = 2 * nnode;
    for i in 0..nnode {
        let c = 2 * i;
        let dx = dndx[i][0];
        let dy = dndx[i][1];
        b[0 * n + c] = dx;
        b[1 * n + c + 1] = dy;
        b[2 * n + c] = dy;
        b[2 * n + c + 1] = dx;
    }
}

const HEX_XI: [[f64; 3]; 8] = [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
];

fn hex8_shape(xi: f64, eta: f64, zeta: f64) -> ([f64; 8], [[f64; 3]; 8]) {
    let mut n = [0.0; 8];
    let mut dn = [[0.0; 3]; 8];
    for i in 0..8 {
        let x = HEX_XI[i][0];
        let e = HEX_XI[i][1];
        let z = HEX_XI[i][2];
        n[i] = 0.125 * (1.0 + x * xi) * (1.0 + e * eta) * (1.0 + z * zeta);
        dn[i][0] = 0.125 * x * (1.0 + e * eta) * (1.0 + z * zeta);
        dn[i][1] = 0.125 * e * (1.0 + x * xi) * (1.0 + z * zeta);
        dn[i][2] = 0.125 * z * (1.0 + x * xi) * (1.0 + e * eta);
    }
    (n, dn)
}

fn hex8_dndx(xyz: &[[f64; 3]; 8], xi: f64, eta: f64, zeta: f64) -> Result<([[f64; 3]; 8], f64, [f64; 8])> {
    let (n, dn) = hex8_shape(xi, eta, zeta);
    let mut j = [[0.0; 3]; 3];
    for a in 0..8 {
        for p in 0..3 {
            for q in 0..3 {
                j[q][p] += dn[a][p] * xyz[a][q];
            }
        }
    }
    // J_{q p} = d x_q / d ξ_p
    let (inv, det) = invert3(j)?;
    let mut dndx = [[0.0; 3]; 8];
    for a in 0..8 {
        for i in 0..3 {
            dndx[a][i] = inv[0][i] * dn[a][0] + inv[1][i] * dn[a][1] + inv[2][i] * dn[a][2];
        }
    }
    Ok((dndx, det, n))
}

fn hex8_stiffness(xyz: &[[f64; 3]; 8], e: f64, nu: f64) -> Result<(Vec<f64>, f64)> {
    let d = d_iso_3d(e, nu)?;
    let n = 24usize;
    let mut ke = vec![0.0; n * n];
    let mut vol = 0.0;
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            for &zeta in &pts {
                let (dndx, det, _) = hex8_dndx(xyz, xi, eta, zeta)?;
                if det <= 0.0 {
                    return err("C3D8: negative Jakobideterminante.");
                }
                let mut b = vec![0.0; 6 * n];
                fill_b3(&mut b, 8, &dndx);
                gemm_bt_d_b(&mut ke, n, &b, 6, &d, det);
                vol += det;
            }
        }
    }
    Ok((ke, vol))
}

fn tet4_stiffness(xyz: &[[f64; 3]; 4], e: f64, nu: f64) -> Result<(Vec<f64>, f64, [[f64; 3]; 4])> {
    let mut j = [[0.0; 3]; 3];
    for p in 0..3 {
        for q in 0..3 {
            j[q][p] = xyz[p + 1][q] - xyz[0][q];
        }
    }
    let (inv, det) = invert3(j)?;
    if det <= 0.0 {
        return err("C3D4: negative Jakobideterminante (Knotenreihenfolge).");
    }
    let vol = det / 6.0;
    // dN/dξ for parent tet N0=1-r-s-t, N1=r, N2=s, N3=t
    let dn = [[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut dndx = [[0.0; 3]; 4];
    for a in 0..4 {
        for i in 0..3 {
            dndx[a][i] = inv[0][i] * dn[a][0] + inv[1][i] * dn[a][1] + inv[2][i] * dn[a][2];
        }
    }
    let d = d_iso_3d(e, nu)?;
    let n = 12usize;
    let mut ke = vec![0.0; n * n];
    let mut b = vec![0.0; 6 * n];
    fill_b3(&mut b, 4, &dndx);
    gemm_bt_d_b(&mut ke, n, &b, 6, &d, vol);
    Ok((ke, vol, dndx))
}

const QUAD_XI: [[f64; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];

fn quad4_shape(xi: f64, eta: f64) -> ([f64; 4], [[f64; 2]; 4]) {
    let mut n = [0.0; 4];
    let mut dn = [[0.0; 2]; 4];
    for i in 0..4 {
        let x = QUAD_XI[i][0];
        let e = QUAD_XI[i][1];
        n[i] = 0.25 * (1.0 + x * xi) * (1.0 + e * eta);
        dn[i][0] = 0.25 * x * (1.0 + e * eta);
        dn[i][1] = 0.25 * e * (1.0 + x * xi);
    }
    (n, dn)
}

fn quad4_dndx(xy: &[[f64; 2]; 4], xi: f64, eta: f64) -> Result<([[f64; 2]; 4], f64, [f64; 4])> {
    let (n, dn) = quad4_shape(xi, eta);
    let mut j = [[0.0; 2]; 2];
    for a in 0..4 {
        j[0][0] += dn[a][0] * xy[a][0];
        j[0][1] += dn[a][1] * xy[a][0];
        j[1][0] += dn[a][0] * xy[a][1];
        j[1][1] += dn[a][1] * xy[a][1];
    }
    let (inv, det) = invert2(j)?;
    let mut dndx = [[0.0; 2]; 4];
    for a in 0..4 {
        dndx[a][0] = inv[0][0] * dn[a][0] + inv[1][0] * dn[a][1];
        dndx[a][1] = inv[0][1] * dn[a][0] + inv[1][1] * dn[a][1];
    }
    Ok((dndx, det, n))
}

fn quad4_stiffness(xy: &[[f64; 2]; 4], e: f64, nu: f64, t: f64, plane_strain: bool) -> Result<(Vec<f64>, f64)> {
    let d = if plane_strain {
        d_plane_strain(e, nu)?
    } else {
        d_plane_stress(e, nu)?
    };
    let n = 8usize;
    let mut ke = vec![0.0; n * n];
    let mut area = 0.0;
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            let (dndx, det, _) = quad4_dndx(xy, xi, eta)?;
            if det <= 0.0 {
                return err("CPS4/CPE4: negative Jakobideterminante.");
            }
            let mut b = vec![0.0; 3 * n];
            fill_b2(&mut b, 4, &dndx);
            gemm_bt_d_b(&mut ke, n, &b, 3, &d, t * det);
            area += det;
        }
    }
    Ok((ke, area))
}

fn tri3_stiffness(xy: &[[f64; 2]; 3], e: f64, nu: f64, t: f64, plane_strain: bool) -> Result<(Vec<f64>, f64, [[f64; 2]; 3])> {
    let x1 = xy[0][0];
    let y1 = xy[0][1];
    let x2 = xy[1][0];
    let y2 = xy[1][1];
    let x3 = xy[2][0];
    let y3 = xy[2][1];
    let two_a = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
    if two_a <= 0.0 {
        return err("CPS3/CPE3: nicht-positive Fläche (Knotenreihenfolge).");
    }
    let a = 0.5 * two_a;
    let mut dndx = [[0.0; 2]; 3];
    dndx[0] = [(y2 - y3) / two_a, (x3 - x2) / two_a];
    dndx[1] = [(y3 - y1) / two_a, (x1 - x3) / two_a];
    dndx[2] = [(y1 - y2) / two_a, (x2 - x1) / two_a];
    let d = if plane_strain {
        d_plane_strain(e, nu)?
    } else {
        d_plane_stress(e, nu)?
    };
    let n = 6usize;
    let mut ke = vec![0.0; n * n];
    let mut b = vec![0.0; 3 * n];
    fill_b2(&mut b, 3, &dndx);
    gemm_bt_d_b(&mut ke, n, &b, 3, &d, t * a);
    Ok((ke, a, dndx))
}

pub struct KeFe {
    pub ke: Vec<f64>,
    pub fe: Vec<f64>,
    pub ndof: usize,
    pub volume: f64,
}

pub fn element_ke(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    e: f64,
    nu: f64,
    thickness: f64,
    section: Option<&BeamSection>,
) -> Result<KeFe> {
    match kind {
        ElemKind::Hex8 => {
            let mut p = [[0.0; 3]; 8];
            for i in 0..8 {
                p[i] = xyz[i];
            }
            let (ke, vol) = hex8_stiffness(&p, e, nu)?;
            let ndof = 24;
            Ok(KeFe {
                ke,
                fe: vec![0.0; ndof],
                ndof,
                volume: vol,
            })
        }
        ElemKind::Tet4 => {
            let mut p = [[0.0; 3]; 4];
            for i in 0..4 {
                p[i] = xyz[i];
            }
            let (ke, vol, _) = tet4_stiffness(&p, e, nu)?;
            Ok(KeFe {
                ke,
                fe: vec![0.0; 12],
                ndof: 12,
                volume: vol,
            })
        }
        ElemKind::Quad4Ps | ElemKind::Quad4Pe => {
            let mut p = [[0.0; 2]; 4];
            for i in 0..4 {
                p[i] = [xyz[i][0], xyz[i][1]];
            }
            let (ke, area) = quad4_stiffness(&p, e, nu, thickness, kind.is_plane_strain())?;
            Ok(KeFe {
                ke,
                fe: vec![0.0; 8],
                ndof: 8,
                volume: area * thickness,
            })
        }
        ElemKind::Tri3Ps | ElemKind::Tri3Pe => {
            let mut p = [[0.0; 2]; 3];
            for i in 0..3 {
                p[i] = [xyz[i][0], xyz[i][1]];
            }
            let (ke, area, _) = tri3_stiffness(&p, e, nu, thickness, kind.is_plane_strain())?;
            Ok(KeFe {
                ke,
                fe: vec![0.0; 6],
                ndof: 6,
                volume: area * thickness,
            })
        }
        ElemKind::Beam31 | ElemKind::Beam32 => {
            let sec = section.ok_or_else(|| {
                crate::error::FemError("Balkenelement ohne *BEAM SECTION.".into())
            })?;
            let (ke, len) = beam::stiffness(kind, xyz, e, nu, sec)?;
            let ndof = 6 * kind.nnodes();
            Ok(KeFe {
                ke,
                fe: vec![0.0; ndof],
                ndof,
                volume: len * sec.area,
            })
        }
        ElemKind::Hex20 | ElemKind::Hex20R => {
            let (ke, vol) = quadratic::hex20_stiffness(xyz, e, nu, kind.reduced_int())?;
            Ok(KeFe {
                ke,
                fe: vec![0.0; 60],
                ndof: 60,
                volume: vol,
            })
        }
        ElemKind::Tet10 => {
            let (ke, vol) = quadratic::tet10_stiffness(xyz, e, nu)?;
            Ok(KeFe {
                ke,
                fe: vec![0.0; 30],
                ndof: 30,
                volume: vol,
            })
        }
        ElemKind::Quad8Ps | ElemKind::Quad8Pe | ElemKind::Quad8RPs | ElemKind::Quad8RPe => {
            let mut p = [[0.0; 2]; 8];
            for i in 0..8 {
                p[i] = [xyz[i][0], xyz[i][1]];
            }
            let (ke, area) = quadratic::quad8_stiffness(
                &p,
                e,
                nu,
                thickness,
                kind.is_plane_strain(),
                kind.reduced_int(),
            )?;
            Ok(KeFe {
                ke,
                fe: vec![0.0; 16],
                ndof: 16,
                volume: area * thickness,
            })
        }
        ElemKind::Tri6Ps | ElemKind::Tri6Pe => {
            let mut p = [[0.0; 2]; 6];
            for i in 0..6 {
                p[i] = [xyz[i][0], xyz[i][1]];
            }
            let (ke, area) =
                quadratic::tri6_stiffness(&p, e, nu, thickness, kind.is_plane_strain())?;
            Ok(KeFe {
                ke,
                fe: vec![0.0; 12],
                ndof: 12,
                volume: area * thickness,
            })
        }
        ElemKind::Shell4
        | ElemKind::Shell4R
        | ElemKind::Shell3
        | ElemKind::Shell8
        | ElemKind::Shell8R
        | ElemKind::Shell6 => {
            let (ke, area) = shell::stiffness(kind, xyz, e, nu, thickness)?;
            let ndof = 6 * kind.nnodes();
            Ok(KeFe {
                ke,
                fe: vec![0.0; ndof],
                ndof,
                volume: area * thickness,
            })
        }
    }
}

/// Stress at nodes of an element from displacement ue (element dof vector).
/// Returns Voigt stress [sxx,syy,szz,sxy,syz,szx] per element node (extrapolated / constant).
pub fn element_nodal_stress(
    kind: ElemKind,
    xyz: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    section: Option<&BeamSection>,
    thickness: f64,
) -> Result<Vec<[f64; 6]>> {
    match kind {
        ElemKind::Hex8 => hex8_nodal_stress(xyz, ue, e, nu),
        ElemKind::Tet4 => tet4_nodal_stress(xyz, ue, e, nu),
        ElemKind::Quad4Ps | ElemKind::Quad4Pe => {
            quad4_nodal_stress(xyz, ue, e, nu, kind.is_plane_strain())
        }
        ElemKind::Tri3Ps | ElemKind::Tri3Pe => {
            tri3_nodal_stress(xyz, ue, e, nu, kind.is_plane_strain())
        }
        ElemKind::Beam31 | ElemKind::Beam32 => {
            let sec = section.ok_or_else(|| {
                crate::error::FemError("Balkenelement ohne *BEAM SECTION.".into())
            })?;
            beam::nodal_stress(kind, xyz, ue, e, nu, sec)
        }
        ElemKind::Hex20 | ElemKind::Hex20R => quadratic::hex20_nodal_stress(xyz, ue, e, nu),
        ElemKind::Tet10 => quadratic::tet10_nodal_stress(xyz, ue, e, nu),
        ElemKind::Quad8Ps | ElemKind::Quad8Pe | ElemKind::Quad8RPs | ElemKind::Quad8RPe => {
            quadratic::quad8_nodal_stress(xyz, ue, e, nu, kind.is_plane_strain())
        }
        ElemKind::Tri6Ps | ElemKind::Tri6Pe => {
            quadratic::tri6_nodal_stress(xyz, ue, e, nu, kind.is_plane_strain())
        }
        ElemKind::Shell4
        | ElemKind::Shell4R
        | ElemKind::Shell3
        | ElemKind::Shell8
        | ElemKind::Shell8R
        | ElemKind::Shell6 => shell::nodal_stress(kind, xyz, ue, e, nu, thickness),
    }
}

pub(crate) fn sigma_from_b(b: &[f64], nrow: usize, n: usize, d: &[f64], ue: &[f64]) -> Vec<f64> {
    let mut eps = vec![0.0; nrow];
    for r in 0..nrow {
        let mut s = 0.0;
        for j in 0..n {
            s += b[r * n + j] * ue[j];
        }
        eps[r] = s;
    }
    let mut sig = vec![0.0; nrow];
    for r in 0..nrow {
        let mut s = 0.0;
        for k in 0..nrow {
            s += d[r * nrow + k] * eps[k];
        }
        sig[r] = s;
    }
    sig
}

fn hex8_nodal_stress(xyz: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<Vec<[f64; 6]>> {
    let mut p = [[0.0; 3]; 8];
    for i in 0..8 {
        p[i] = xyz[i];
    }
    let d = d_iso_3d(e, nu)?;
    let n = 24usize;
    let mut gsig = Vec::with_capacity(8);
    for i in 0..8 {
        let xi = HEX_XI[i][0] * G2;
        let eta = HEX_XI[i][1] * G2;
        let zeta = HEX_XI[i][2] * G2;
        let (dndx, _, _) = hex8_dndx(&p, xi, eta, zeta)?;
        let mut b = vec![0.0; 6 * n];
        fill_b3(&mut b, 8, &dndx);
        let s = sigma_from_b(&b, 6, n, &d, ue);
        let mut six = [0.0; 6];
        six.copy_from_slice(&s);
        gsig.push(six);
    }
    let s3 = 3.0_f64.sqrt();
    let mut out = vec![[0.0; 6]; 8];
    for a in 0..8 {
        let xi = HEX_XI[a][0] * s3;
        let eta = HEX_XI[a][1] * s3;
        let zeta = HEX_XI[a][2] * s3;
        let (nshp, _) = hex8_shape(xi, eta, zeta);
        for c in 0..6 {
            let mut v = 0.0;
            for g in 0..8 {
                v += nshp[g] * gsig[g][c];
            }
            out[a][c] = v;
        }
    }
    Ok(out)
}

fn tet4_nodal_stress(xyz: &[[f64; 3]], ue: &[f64], e: f64, nu: f64) -> Result<Vec<[f64; 6]>> {
    let mut p = [[0.0; 3]; 4];
    for i in 0..4 {
        p[i] = xyz[i];
    }
    let (_, _, dndx) = tet4_stiffness(&p, e, nu)?;
    let d = d_iso_3d(e, nu)?;
    let n = 12usize;
    let mut b = vec![0.0; 6 * n];
    fill_b3(&mut b, 4, &dndx);
    let s = sigma_from_b(&b, 6, n, &d, ue);
    let mut six = [0.0; 6];
    six.copy_from_slice(&s);
    Ok(vec![six; 4])
}

fn to6_from_plane(s: &[f64], e: f64, nu: f64, plane_strain: bool) -> [f64; 6] {
    let sxx = s[0];
    let syy = s[1];
    let sxy = s[2];
    let szz = if plane_strain {
        nu * (sxx + syy)
    } else {
        0.0
    };
    let _ = e;
    [sxx, syy, szz, sxy, 0.0, 0.0]
}

fn quad4_nodal_stress(
    xyz: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    plane_strain: bool,
) -> Result<Vec<[f64; 6]>> {
    let mut p = [[0.0; 2]; 4];
    for i in 0..4 {
        p[i] = [xyz[i][0], xyz[i][1]];
    }
    let d = if plane_strain {
        d_plane_strain(e, nu)?
    } else {
        d_plane_stress(e, nu)?
    };
    let n = 8usize;
    let mut gsig = Vec::with_capacity(4);
    for i in 0..4 {
        let xi = QUAD_XI[i][0] * G2;
        let eta = QUAD_XI[i][1] * G2;
        let (dndx, _, _) = quad4_dndx(&p, xi, eta)?;
        let mut b = vec![0.0; 3 * n];
        fill_b2(&mut b, 4, &dndx);
        let s = sigma_from_b(&b, 3, n, &d, ue);
        gsig.push(s);
    }
    let s3 = 3.0_f64.sqrt();
    let mut out = vec![[0.0; 6]; 4];
    for a in 0..4 {
        let xi = QUAD_XI[a][0] * s3;
        let eta = QUAD_XI[a][1] * s3;
        let (nshp, _) = quad4_shape(xi, eta);
        let mut sp = [0.0; 3];
        for c in 0..3 {
            for g in 0..4 {
                sp[c] += nshp[g] * gsig[g][c];
            }
        }
        out[a] = to6_from_plane(&sp, e, nu, plane_strain);
    }
    Ok(out)
}

fn tri3_nodal_stress(
    xyz: &[[f64; 3]],
    ue: &[f64],
    e: f64,
    nu: f64,
    plane_strain: bool,
) -> Result<Vec<[f64; 6]>> {
    let mut p = [[0.0; 2]; 3];
    for i in 0..3 {
        p[i] = [xyz[i][0], xyz[i][1]];
    }
    let (_, _, dndx) = tri3_stiffness(&p, e, nu, 1.0, plane_strain)?;
    let d = if plane_strain {
        d_plane_strain(e, nu)?
    } else {
        d_plane_stress(e, nu)?
    };
    let n = 6usize;
    let mut b = vec![0.0; 3 * n];
    fill_b2(&mut b, 3, &dndx);
    let s = sigma_from_b(&b, 3, n, &d, ue);
    let six = to6_from_plane(&s, e, nu, plane_strain);
    Ok(vec![six; 3])
}

pub fn von_mises(s: &[f64; 6]) -> f64 {
    let sx = s[0];
    let sy = s[1];
    let sz = s[2];
    let txy = s[3];
    let tyz = s[4];
    let tzx = s[5];
    (0.5 * ((sx - sy).powi(2) + (sy - sz).powi(2) + (sz - sx).powi(2))
        + 3.0 * (txy * txy + tyz * tyz + tzx * tzx))
        .sqrt()
}

/// Face pressure on C3D8. face 1..=6 (CalculiX P1..P6). Positive pressure into the element.
pub fn hex8_face_pressure(xyz: &[[f64; 3]; 8], face: i32, p: f64) -> Result<[f64; 24]> {
    // Local 0-based node indices per Abaqus/CalculiX face
    let faces: [[usize; 4]; 6] = [
        [0, 1, 2, 3],
        [4, 7, 6, 5],
        [0, 4, 5, 1],
        [1, 5, 6, 2],
        [2, 6, 7, 3],
        [3, 7, 4, 0],
    ];
    if !(1..=6).contains(&face) {
        return err(format!("Ungültige C3D8-Fläche P{face}"));
    }
    let fi = (face - 1) as usize;
    let mut fe = [0.0; 24];
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            let (n, dnxi, dneta) = {
                let mut n = [0.0; 4];
                let mut dnxi = [0.0; 4];
                let mut dneta = [0.0; 4];
                for i in 0..4 {
                    let x = QUAD_XI[i][0];
                    let e = QUAD_XI[i][1];
                    n[i] = 0.25 * (1.0 + x * xi) * (1.0 + e * eta);
                    dnxi[i] = 0.25 * x * (1.0 + e * eta);
                    dneta[i] = 0.25 * e * (1.0 + x * xi);
                }
                (n, dnxi, dneta)
            };
            let mut rxi = [0.0; 3];
            let mut reta = [0.0; 3];
            for a in 0..4 {
                let q = xyz[faces[fi][a]];
                for k in 0..3 {
                    rxi[k] += dnxi[a] * q[k];
                    reta[k] += dneta[a] * q[k];
                }
            }
            let nx = rxi[1] * reta[2] - rxi[2] * reta[1];
            let ny = rxi[2] * reta[0] - rxi[0] * reta[2];
            let nz = rxi[0] * reta[1] - rxi[1] * reta[0];
            // traction = -p * n_outward, dA vector is n_outward * dξ dη
            let tx = -p * nx;
            let ty = -p * ny;
            let tz = -p * nz;
            for a in 0..4 {
                let gd = 3 * faces[fi][a];
                fe[gd] += n[a] * tx;
                fe[gd + 1] += n[a] * ty;
                fe[gd + 2] += n[a] * tz;
            }
        }
    }
    Ok(fe)
}

/// Edge pressure on 2D quad, P1..P4, positive into the element.
pub fn quad4_edge_pressure(xy: &[[f64; 2]; 4], face: i32, p: f64, thickness: f64) -> Result<[f64; 8]> {
    let edges: [[usize; 2]; 4] = [[0, 1], [1, 2], [2, 3], [3, 0]];
    if !(1..=4).contains(&face) {
        return err(format!("Ungültige CPS4-Kante P{face}"));
    }
    let e = (face - 1) as usize;
    let a = edges[e][0];
    let b = edges[e][1];
    let dx = xy[b][0] - xy[a][0];
    let dy = xy[b][1] - xy[a][1];
    // inward normal for CCW: (-dy, dx) wait rotate tangent 90 CCW = (-dy, dx)
    // tangent (dx,dy), inward (CCW interior is left): (-dy, dx)
    let nx = -dy;
    let ny = dx;
    // traction = p * n_inward; n_inward dL = (-dy, dx) which already includes length
    // each node gets half
    let fx = p * thickness * nx * 0.5;
    let fy = p * thickness * ny * 0.5;
    let mut fe = [0.0; 8];
    fe[2 * a] += fx;
    fe[2 * a + 1] += fy;
    fe[2 * b] += fx;
    fe[2 * b + 1] += fy;
    Ok(fe)
}

pub fn hex8_body_force(xyz: &[[f64; 3]; 8], bx: f64, by: f64, bz: f64) -> Result<[f64; 24]> {
    let mut fe = [0.0; 24];
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            for &zeta in &pts {
                let (_, det, n) = hex8_dndx(xyz, xi, eta, zeta)?;
                for a in 0..8 {
                    fe[3 * a] += n[a] * bx * det;
                    fe[3 * a + 1] += n[a] * by * det;
                    fe[3 * a + 2] += n[a] * bz * det;
                }
            }
        }
    }
    Ok(fe)
}

pub fn tet4_body_force(xyz: &[[f64; 3]; 4], bx: f64, by: f64, bz: f64) -> Result<[f64; 12]> {
    let mut j = [[0.0; 3]; 3];
    for p in 0..3 {
        for q in 0..3 {
            j[q][p] = xyz[p + 1][q] - xyz[0][q];
        }
    }
    let (_, det) = invert3(j)?;
    let vol = det / 6.0;
    let mut fe = [0.0; 12];
    for a in 0..4 {
        fe[3 * a] = vol * 0.25 * bx;
        fe[3 * a + 1] = vol * 0.25 * by;
        fe[3 * a + 2] = vol * 0.25 * bz;
    }
    Ok(fe)
}

pub fn quad4_body_force(
    xy: &[[f64; 2]; 4],
    bx: f64,
    by: f64,
    thickness: f64,
) -> Result<[f64; 8]> {
    let mut fe = [0.0; 8];
    let pts = [-G2, G2];
    for &xi in &pts {
        for &eta in &pts {
            let (_, det, n) = quad4_dndx(xy, xi, eta)?;
            for a in 0..4 {
                fe[2 * a] += n[a] * bx * det * thickness;
                fe[2 * a + 1] += n[a] * by * det * thickness;
            }
        }
    }
    Ok(fe)
}

pub fn xyz_of(kind: ElemKind, all: &[[f64; 3]], conn: &[usize]) -> Vec<[f64; 3]> {
    let n = kind.nnodes();
    (0..n).map(|i| all[conn[i]]).collect()
}
