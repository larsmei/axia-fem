//! Constitutive helpers used by Newton (1.3) and the T3D2 1D return map.

/// Piecewise-linear yield curve: CalculiX `*PLASTIC` is (sy, peeq).
/// Stored as (peeq, sy).
pub fn yield_from_curve(curve: &[(f64, f64)], peeq: f64) -> (f64, f64) {
    if curve.is_empty() {
        return (0.0, 0.0);
    }
    if peeq <= curve[0].0 {
        return (curve[0].1, 0.0);
    }
    for w in curve.windows(2) {
        let (p0, s0) = w[0];
        let (p1, s1) = w[1];
        if peeq <= p1 {
            let t = if (p1 - p0).abs() < 1e-16 {
                0.0
            } else {
                (peeq - p0) / (p1 - p0)
            };
            let sy = s0 + t * (s1 - s0);
            let h = if (p1 - p0).abs() < 1e-16 {
                0.0
            } else {
                (s1 - s0) / (p1 - p0)
            };
            return (sy, h.max(0.0));
        }
    }
    let last = curve[curve.len() - 1];
    (last.1, 0.0)
}

/// 1D J2 with isotropic hardening. Returns (sigma, tangent E_t).
pub fn truss_1d_stress(e: f64, eps: f64, curve: Option<&[(f64, f64)]>) -> (f64, f64) {
    let Some(c) = curve else {
        return (e * eps, e);
    };
    if c.is_empty() {
        return (e * eps, e);
    }
    let sign = if eps >= 0.0 { 1.0 } else { -1.0 };
    let aeps = eps.abs();
    let mut pe = 0.0;
    let (sy0, _) = yield_from_curve(c, 0.0);
    if aeps <= sy0 / e {
        return (sign * e * aeps, e);
    }
    for _ in 0..40 {
        let (sy, h) = yield_from_curve(c, pe);
        let eps_of = sy / e + pe;
        let r = aeps - eps_of;
        if r.abs() < 1e-14 {
            return (sign * sy, (e * h) / (e + h).max(1e-30));
        }
        let den = h / e + 1.0;
        pe += r / den.max(1e-30);
        pe = pe.max(0.0);
    }
    let (sy, h) = yield_from_curve(c, pe);
    (sign * sy, (e * h) / (e + h).max(1e-30))
}

pub struct NewtonCtrl {
    pub max_iter: usize,
    pub rtol: f64,
}

impl Default for NewtonCtrl {
    fn default() -> Self {
        Self {
            max_iter: 25,
            rtol: 1e-8,
        }
    }
}
