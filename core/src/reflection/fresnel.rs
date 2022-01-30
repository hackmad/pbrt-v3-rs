//! Fresnel Dielectrics and Conductors

#![allow(dead_code)]

use super::*;
use bumpalo::Bump;
use std::mem::swap;

/// Interface for computing Fresnel reflection coefficients.
pub enum Fresnel<'arena> {
    NoOp(&'arena mut FresnelNoOp),
    Dielectric(&'arena mut FresnelDielectric),
    Conductor(&'arena mut FresnelConductor),
}

impl<'arena> Fresnel<'arena> {
    /// Clone into a newly allocated instance.
    ///
    /// * `arena` - The memory arena used for allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc(&self, arena: &'arena Bump) -> &'arena mut Fresnel<'arena> {
        match self {
            Self::NoOp(f) => f.clone_alloc(arena),
            Self::Dielectric(f) => f.clone_alloc(arena),
            Self::Conductor(f) => f.clone_alloc(arena),
        }
    }

    /// Returns the amount of light reflected by the surface.
    ///
    /// * `cos_thata_i` - Cosine of the angle made by incident direction and
    ///                   surface normal.
    pub fn evaluate(&self, cos_theta_i: Float) -> Spectrum {
        match self {
            Self::NoOp(f) => f.evaluate(cos_theta_i),
            Self::Dielectric(f) => f.evaluate(cos_theta_i),
            Self::Conductor(f) => f.evaluate(cos_theta_i),
        }
    }
}

/// Implements `Fresnel` for dielectric materials.
#[derive(Default)]
pub struct FresnelDielectric {
    /// Index of refraction for exterior side of the surface.
    eta_i: Float,

    /// Index of refraction for interior side of the surface.
    eta_t: Float,
}

impl FresnelDielectric {
    /// Allocate a new `FresnelDielectric`.
    ///
    /// * `arena` - The arena for memory allocations.
    /// * `eta_i` - Index of refraction for exterior side of the surface.
    /// * `eta_t` - Index of refraction for interior side of the surface.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc<'arena>(arena: &'arena Bump, eta_i: Float, eta_t: Float) -> &'arena mut Fresnel {
        let f = arena.alloc(Self { eta_i, eta_t });
        arena.alloc(Fresnel::Dielectric(f))
    }

    /// Clone into a newly allocated a new instance of `FresnelDielectric`.
    ///
    /// * `arena` - The arena for memory allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc<'arena>(&self, arena: &'arena Bump) -> &'arena mut Fresnel<'arena> {
        let f = arena.alloc(Self {
            eta_i: self.eta_i,
            eta_t: self.eta_t,
        });
        arena.alloc(Fresnel::Dielectric(f))
    }

    /// Returns the amount of light reflected by the surface.
    ///
    /// * `cos_thata_i` - Cosine of the angle made by incident direction and
    ///                   surface normal.
    pub fn evaluate(&self, cos_theta_i: Float) -> Spectrum {
        Spectrum::new(fr_dielectric(cos_theta_i, self.eta_i, self.eta_t))
    }
}

/// Implements `Fresnel` for conductors materials.
#[derive(Default)]
pub struct FresnelConductor {
    /// Index of refraction for exterior side of the surface.
    eta_i: Spectrum,

    /// Index of refraction for interior side of the surface.
    eta_t: Spectrum,

    /// Absorption coefficient.
    k: Spectrum,
}

impl FresnelConductor {
    /// Allocate a new `FresnelConductor`.
    ///
    /// * `arena` - The arena for memory allocations.
    /// * `eta_i` - Index of refraction for exterior side of the surface.
    /// * `eta_t` - Index of refraction for interior side of the surface.
    /// * `k`     - Absorption coefficient.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc<'arena>(
        arena: &'arena Bump,
        eta_i: Spectrum,
        eta_t: Spectrum,
        k: Spectrum,
    ) -> &'arena mut Fresnel {
        let f = arena.alloc(Self { eta_i, eta_t, k });
        arena.alloc(Fresnel::Conductor(f))
    }

    /// Clone into a newly allocated a new instance of `FresnelConductor`.
    ///
    /// * `arena` - The arena for memory allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc<'arena>(&self, arena: &'arena Bump) -> &'arena mut Fresnel<'arena> {
        let f = arena.alloc(Self {
            eta_i: self.eta_i,
            eta_t: self.eta_t,
            k: self.k,
        });
        arena.alloc(Fresnel::Conductor(f))
    }

    /// Returns the amount of light reflected by the surface.
    ///
    /// * `cos_thata_i` - Cosine of the angle made by incident direction and
    ///                   surface normal.
    pub fn evaluate(&self, cos_theta_i: Float) -> Spectrum {
        // Need to take abs(cos_theta_i)) so angle is measured on same side as
        // normal.
        fr_conductor(abs(cos_theta_i), self.eta_i, self.eta_t, self.k)
    }
}

/// Implements `Fresnel` for materials that reflect 100% of all incoming light.
#[derive(Default)]
pub struct FresnelNoOp {}
impl FresnelNoOp {
    /// Allocator a new `FresnelNoOp`.
    ///
    /// * `allocator` - The allocator.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc<'arena>(arena: &'arena Bump) -> &'arena mut Fresnel {
        let f = arena.alloc(Self {});
        arena.alloc(Fresnel::NoOp(f))
    }

    /// Clone into a newly allocated a new instance of `FresnelNoOp`.
    ///
    /// * `arena` - The arena for memory allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc<'arena>(&self, arena: &'arena Bump) -> &'arena mut Fresnel<'arena> {
        let f = arena.alloc(Self {});
        arena.alloc(Fresnel::NoOp(f))
    }

    /// Returns the amount of light reflected by the surface.
    ///
    /// * `cos_thata_i` - Cosine of the angle made by incident direction and
    ///                   surface normal.
    pub fn evaluate(&self, _cos_theta_i: Float) -> Spectrum {
        Spectrum::ONE
    }
}

/// Returns the fresnel reflection for dielectric materials and unpolarized light.
///
/// * `cos_theta_i` - cos(θi) for angle between incident direction and geometric
///                   surface normal.
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
    let sin_theta_i = max(0.0, 1.0 - cos_theta_i * cos_theta_i).sqrt();
    let sin_theta_t = eta_i / eta_t * sin_theta_i;

    // Handle total internal reflection.
    if sin_theta_t >= 1.0 {
        1.0
    } else {
        let cos_theta_t = max(0.0, 1.0 - sin_theta_t * sin_theta_t).sqrt();
        let r_parl = ((eta_t * cos_theta_i) - (eta_i * cos_theta_t))
            / ((eta_t * cos_theta_i) + (eta_i * cos_theta_t));
        let r_perp = ((eta_i * cos_theta_i) - (eta_t * cos_theta_t))
            / ((eta_i * cos_theta_i) + (eta_t * cos_theta_t));
        (r_parl * r_parl + r_perp * r_perp) / 2.0
    }
}

/// Returns the Fresnel reflection at the boundary between a conductor and
/// dielectric medium for unpolarized light.
///
/// * `cos_theta_i` - cos(θi) for angle between incident direction and geometric
///                   surface normal on the same side as incident direction `wi`.
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
