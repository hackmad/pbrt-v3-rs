//! Specular Transmission

#![allow(dead_code)]

use super::*;
use crate::core::material::*;

/// BTDF for physically plausible specular transmission using Fresnel interface.
#[derive(Copy, Clone)]
pub struct SpecularTransmission {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Fresnel interface for dielectrics and conductors.
    fresnel: FresnelDielectric,

    /// Spectrum used to scale the transmitted colour.
    t: Spectrum,

    /// Index of refraction above the surface (same side as surface normal).
    eta_a: Float,

    /// Index of refraction below the surface (opposite side as surface normal).
    eta_b: Float,

    /// Indicates whether incident ray started from a light source or from camera.
    mode: TransportMode,
}

impl SpecularTransmission {
    /// Create a new instance of `SpecularTransmission`.
    ///
    /// * `t`     - Spectrum used to scale the transmitted colour.
    /// * `eta_a` - Index of refraction above the surface (same side as surface
    ///             normal).
    /// * `eta_b` - Index of refraction below the surface (opposite side as surface
    ///             normal).
    /// * `mode`  - Indicates whether incident ray started from a light source
    ///             or from camera.
    pub fn new(t: Spectrum, eta_a: Float, eta_b: Float, mode: TransportMode) -> Self {
        Self {
            bxdf_type: BxDFType::from(BSDF_TRANSMISSION | BSDF_SPECULAR),
            fresnel: FresnelDielectric::new(eta_a, eta_b),
            t,
            eta_a,
            eta_b,
            mode,
        }
    }
}

impl BxDF for SpecularTransmission {
    /// Returns the BxDF type.
    fn get_type(&self) -> BxDFType {
        self.bxdf_type
    }

    /// Returns the value of the distribution function for the given pair of
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        // No scattering is returned.
        Spectrum::new(0.0)
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    fn sample_f(&self, wo: &Vector3f, _u: &Point2f) -> (Spectrum, Float, Vector3f, BxDFType) {
        // Figure out which $\eta$ is incident and which is transmitted.
        let entering = cos_theta(wo) > 0.0;
        let eta_i = if entering { self.eta_a } else { self.eta_b };
        let eta_t = if entering { self.eta_b } else { self.eta_a };

        // Compute ray direction for specular transmission
        if let Some(wi) = refract(
            wo,
            &Normal3f::new(0.0, 0.0, 1.0).face_forward(wo),
            eta_i / eta_t,
        ) {
            let pdf = 1.0;
            let mut ft = self.t * (Spectrum::new(1.0) - self.fresnel.evaluate(cos_theta(&wi)));

            // Account for non-symmetry with transmission to different medium
            if self.mode == TransportMode::Radiance {
                ft *= (eta_i * eta_i) / (eta_t * eta_t);
            }

            (ft / abs_cos_theta(&wi), pdf, wi, self.bxdf_type)
        } else {
            (Spectrum::new(0.0), 0.0, Vector3f::default(), self.bxdf_type)
        }
    }
}
