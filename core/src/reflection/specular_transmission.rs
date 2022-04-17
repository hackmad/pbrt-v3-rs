//! Specular Transmission

use super::*;
use crate::material::*;
use bumpalo::Bump;
use std::fmt;

/// BTDF for physically plausible specular transmission using Fresnel interface.
pub struct SpecularTransmission<'arena> {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Fresnel interface for dielectrics and conductors.
    fresnel: &'arena mut Fresnel<'arena>,

    /// Spectrum used to scale the transmitted colour.
    t: Spectrum,

    /// Index of refraction above the surface (same side as surface normal).
    eta_a: Float,

    /// Index of refraction below the surface (opposite side as surface normal).
    eta_b: Float,

    /// Indicates whether incident ray started from a light source or from camera.
    mode: TransportMode,
}

impl<'arena> SpecularTransmission<'arena> {
    /// Allocate a new instance of `SpecularTransmission`.
    ///
    /// * `arena`  - The arena for memory allocations.
    /// * `t`      - Spectrum used to scale the transmitted colour.
    /// * `eta_a`  - Index of refraction above the surface (same side as
    ///              surface normal).
    /// * `eta_b`  - Index of refraction below the surface (opposite side as
    ///              surface normal).
    /// * `mode`   - Indicates whether incident ray started from a light
    ///              source or from camera.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc(
        arena: &'arena Bump,
        t: Spectrum,
        eta_a: Float,
        eta_b: Float,
        mode: TransportMode,
    ) -> &'arena mut BxDF<'arena> {
        let fresnel = FresnelDielectric::alloc(arena, eta_a, eta_b);
        let model = arena.alloc(Self {
            bxdf_type: BxDFType::BSDF_TRANSMISSION | BxDFType::BSDF_SPECULAR,
            fresnel,
            t,
            eta_a,
            eta_b,
            mode,
        });
        arena.alloc(BxDF::SpecularTransmission(model))
    }

    /// Clone into a newly allocated a new instance of `SpecularTransmission`.
    ///
    /// * `arena` - The arena for memory allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc(&self, arena: &'arena Bump) -> &'arena mut BxDF<'arena> {
        let fresnel = self.fresnel.clone_alloc(arena);
        let model = arena.alloc(Self {
            bxdf_type: self.bxdf_type,
            fresnel,
            t: self.t,
            eta_a: self.eta_a,
            eta_b: self.eta_b,
            mode: self.mode,
        });
        arena.alloc(BxDF::SpecularTransmission(model))
    }

    /// Returns the BxDF type.
    pub fn get_type(&self) -> BxDFType {
        self.bxdf_type
    }

    /// Returns the value of the distribution function for the given pair of
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        // No scattering is returned.
        Spectrum::ZERO
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    pub fn sample_f(&self, wo: &Vector3f, _u: &Point2f) -> BxDFSample {
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
            let mut ft = self.t * (Spectrum::ONE - self.fresnel.evaluate(cos_theta(&wi)));

            // Account for non-symmetry with transmission to different medium
            if self.mode == TransportMode::Radiance {
                ft *= (eta_i * eta_i) / (eta_t * eta_t);
            }

            BxDFSample::new(ft / abs_cos_theta(&wi), pdf, wi, self.bxdf_type)
        } else {
            BxDFSample::from(self.bxdf_type)
        }
    }

    /// Evaluates the PDF for the sampling method. Default is based on the
    /// cosine-weighted sampling in `BxDF::sample_f()` default implementation.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn pdf(&self, _wo: &Vector3f, _wi: &Vector3f) -> Float {
        0.0
    }
}

impl<'arena> fmt::Display for SpecularTransmission<'arena> {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SpecularTransmission {{ bxdf_type: {}, fresnel: {}, t: {}, eta_a: {}, eta_b: {}, mode: {} }}",
            self.bxdf_type, self.fresnel, self.t, self.eta_a, self.eta_b, self.mode
        )
    }
}
