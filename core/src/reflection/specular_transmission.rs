//! Specular Transmission

use super::*;
use crate::material::*;
use std::fmt;

/// BTDF for physically plausible specular transmission using Fresnel interface.
#[derive(Clone)]
pub struct SpecularTransmission {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Fresnel interface for dielectrics and conductors.
    fresnel: Fresnel,

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
    /// Creates a new instance of `SpecularTransmission`.
    ///
    /// * `t`      - Spectrum used to scale the transmitted colour.
    /// * `eta_a`  - Index of refraction above the surface (same side as surface normal).
    /// * `eta_b`  - Index of refraction below the surface (opposite side as surface normal).
    /// * `mode`   - Indicates whether incident ray started from a light source or from camera.
    pub fn new(t: Spectrum, eta_a: Float, eta_b: Float, mode: TransportMode) -> BxDF {
        let fresnel = FresnelDielectric::new(eta_a, eta_b);
        let model = Self {
            bxdf_type: BxDFType::BSDF_TRANSMISSION | BxDFType::BSDF_SPECULAR,
            fresnel,
            t,
            eta_a,
            eta_b,
            mode,
        };
        BxDF::SpecularTransmission(model)
    }

    /// Returns the BxDF type.
    pub fn get_type(&self) -> BxDFType {
        self.bxdf_type
    }

    /// Returns the value of the distribution function for the given pair of directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        // No scattering is returned.
        Spectrum::ZERO
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    pub fn sample_f(&self, wo: &Vector3f, _u: &Point2f) -> BxDFSample {
        // Figure out which $\eta$ is incident and which is transmitted.
        let entering = cos_theta(wo) > 0.0;
        let eta_i = if entering { self.eta_a } else { self.eta_b };
        let eta_t = if entering { self.eta_b } else { self.eta_a };

        // Compute ray direction for specular transmission
        if let Some(wi) = refract(wo, &Normal3f::new(0.0, 0.0, 1.0).face_forward(wo), eta_i / eta_t) {
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

    /// Evaluates the PDF for the sampling method. Default is based on the cosine-weighted sampling in `BxDF::sample_f()`
    /// default implementation.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn pdf(&self, _wo: &Vector3f, _wi: &Vector3f) -> Float {
        0.0
    }
}

impl fmt::Display for SpecularTransmission {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SpecularTransmission {{ bxdf_type: {}, fresnel: {}, t: {}, eta_a: {}, eta_b: {}, mode: {} }}",
            self.bxdf_type, self.fresnel, self.t, self.eta_a, self.eta_b, self.mode
        )
    }
}
