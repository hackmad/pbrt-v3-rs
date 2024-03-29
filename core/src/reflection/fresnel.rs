//! Fresnel Dielectrics and Conductors

use super::*;
use std::fmt;
use std::mem::swap;

/// Interface for computing Fresnel reflection coefficients.
#[derive(Clone)]
pub enum Fresnel {
    NoOp(FresnelNoOp),
    Dielectric(FresnelDielectric),
    Conductor(FresnelConductor),
}

impl Fresnel {
    /// Returns the amount of light reflected by the surface.
    ///
    /// * `cos_thata_i` - Cosine of the angle made by incident direction and surface normal.
    pub fn evaluate(&self, cos_theta_i: Float) -> Spectrum {
        match self {
            Self::NoOp(f) => f.evaluate(cos_theta_i),
            Self::Dielectric(f) => f.evaluate(cos_theta_i),
            Self::Conductor(f) => f.evaluate(cos_theta_i),
        }
    }
}

impl fmt::Display for Fresnel {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoOp(fr) => write!(f, "Fresnel {{ {} }}", fr),
            Self::Dielectric(fr) => write!(f, "Fresnel {{ {} }}", fr),
            Self::Conductor(fr) => write!(f, "Fresnel {{ {} }}", fr),
        }
    }
}

/// Implements `Fresnel` for dielectric materials.
#[derive(Clone, Default)]
pub struct FresnelDielectric {
    /// Index of refraction for exterior side of the surface.
    eta_i: Float,

    /// Index of refraction for interior side of the surface.
    eta_t: Float,
}

impl FresnelDielectric {
    /// Creates a new `FresnelDielectric`.
    ///
    /// * `eta_i` - Index of refraction for exterior side of the surface.
    /// * `eta_t` - Index of refraction for interior side of the surface.
    pub fn new(eta_i: Float, eta_t: Float) -> Fresnel {
        let f = Self { eta_i, eta_t };
        Fresnel::Dielectric(f)
    }

    /// Returns the amount of light reflected by the surface.
    ///
    /// * `cos_thata_i` - Cosine of the angle made by incident direction and surface normal.
    pub fn evaluate(&self, cos_theta_i: Float) -> Spectrum {
        Spectrum::new(fr_dielectric(cos_theta_i, self.eta_i, self.eta_t))
    }
}

impl fmt::Display for FresnelDielectric {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "FresnelDielectric {{ eta_i: {}, eta_t: {} }}",
            self.eta_i, self.eta_t
        )
    }
}

/// Implements `Fresnel` for conductors materials.
#[derive(Clone, Default)]
pub struct FresnelConductor {
    /// Index of refraction for exterior side of the surface.
    eta_i: Spectrum,

    /// Index of refraction for interior side of the surface.
    eta_t: Spectrum,

    /// Absorption coefficient.
    k: Spectrum,
}

impl FresnelConductor {
    /// Creates a new `FresnelConductor`.
    ///
    /// * `eta_i` - Index of refraction for exterior side of the surface.
    /// * `eta_t` - Index of refraction for interior side of the surface.
    /// * `k`     - Absorption coefficient.
    pub fn new(eta_i: Spectrum, eta_t: Spectrum, k: Spectrum) -> Fresnel {
        let f = Self { eta_i, eta_t, k };
        Fresnel::Conductor(f)
    }

    /// Returns the amount of light reflected by the surface.
    ///
    /// * `cos_thata_i` - Cosine of the angle made by incident direction and surface normal.
    pub fn evaluate(&self, cos_theta_i: Float) -> Spectrum {
        // Need to take abs(cos_theta_i)) so angle is measured on same side as normal.
        fr_conductor(abs(cos_theta_i), self.eta_i, self.eta_t, self.k)
    }
}

impl fmt::Display for FresnelConductor {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "FresnelConductor {{ eta_i: {}, eta_t: {}, k: {} }}",
            self.eta_i, self.eta_t, self.k
        )
    }
}

/// Implements `Fresnel` for materials that reflect 100% of all incoming light.
#[derive(Clone, Default)]
pub struct FresnelNoOp {}
impl FresnelNoOp {
    /// Allocator a new `FresnelNoOp`.
    pub fn new() -> Fresnel {
        let f = Self {};
        Fresnel::NoOp(f)
    }

    /// Returns the amount of light reflected by the surface.
    ///
    /// * `cos_thata_i` - Cosine of the angle made by incident direction and surface normal.
    pub fn evaluate(&self, _cos_theta_i: Float) -> Spectrum {
        Spectrum::ONE
    }
}

impl fmt::Display for FresnelNoOp {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FourierNoOp {{}}")
    }
}

/// Returns the fresnel reflection for dielectric materials and unpolarized light.
///
/// * `cos_theta_i` - cos(θi) for angle between incident direction and geometric surface normal.
/// * `eta_i`       - index of refraction for medium that incident ray is in.
/// * `eta_t`       - index of refraction for medium that incident ray is entering.
pub fn fr_dielectric(cos_theta_i: Float, eta_i: Float, eta_t: Float) -> Float {
    let mut cos_theta_i = clamp(cos_theta_i, -1.0, 1.0);
    let mut eta_i = eta_i;
    let mut eta_t = eta_t;

    // Potentially swap indices of refraction.
    let entering = cos_theta_i > 0.0;
    if !entering {
        swap(&mut eta_i, &mut eta_t);
        cos_theta_i = abs(cos_theta_i);
    }

    // Compute _cosThetaT_ using Snell's law.
    let sin_theta_i = (0.0 as Float).max(1.0 - cos_theta_i * cos_theta_i).sqrt();
    let sin_theta_t = eta_i / eta_t * sin_theta_i;

    // Handle total internal reflection.
    if sin_theta_t >= 1.0 {
        1.0
    } else {
        let cos_theta_t = (0.0 as Float).max(1.0 - sin_theta_t * sin_theta_t).sqrt();
        let r_parl = ((eta_t * cos_theta_i) - (eta_i * cos_theta_t)) / ((eta_t * cos_theta_i) + (eta_i * cos_theta_t));
        let r_perp = ((eta_i * cos_theta_i) - (eta_t * cos_theta_t)) / ((eta_i * cos_theta_i) + (eta_t * cos_theta_t));
        (r_parl * r_parl + r_perp * r_perp) / 2.0
    }
}

/// Returns the Fresnel reflection at the boundary between a conductor and
/// dielectric medium for unpolarized light.
///
/// * `cos_theta_i` - cos(θi) for angle between incident direction and geometric surface normal on the same side as
///                   incident direction `wi`.
/// * `eta_i`       - Index of refraction for medium that incident ray is in.
/// * `eta_t`       - Index of refraction for medium that incident ray is entering.
/// * `k`           - The absorption coefficient.
pub fn fr_conductor(cos_theta_i: Float, eta_i: Spectrum, eta_t: Spectrum, k: Spectrum) -> Spectrum {
    let cos_theta_i = clamp(cos_theta_i, -1.0, 1.0);
    let eta = eta_t / eta_i;
    let eta_k = k / eta_i;

    let cos_theta_i_2 = cos_theta_i * cos_theta_i;
    let sin_theta_i_2 = 1.0 - cos_theta_i;
    let eta_2 = eta * eta;
    let eta_k_2 = eta_k * eta_k;

    let t0 = eta_2 - eta_k_2 - Spectrum::new(sin_theta_i_2);
    let a2_plus_b2 = (t0 * t0 + 4.0 * eta_2 * eta_k_2).sqrt();
    let t1 = a2_plus_b2 + Spectrum::new(cos_theta_i_2);
    let a = (0.5 * (a2_plus_b2 + t0)).sqrt();
    let t2 = 2.0 * cos_theta_i * a;
    let rs = (t1 - t2) / (t1 + t2);

    let t3 = cos_theta_i_2 * a2_plus_b2 + Spectrum::new(sin_theta_i_2 * sin_theta_i_2);
    let t4 = t2 * sin_theta_i_2;
    let rp = rs * (t3 - t4) / (t3 + t4);

    0.5 * (rp + rs)
}
