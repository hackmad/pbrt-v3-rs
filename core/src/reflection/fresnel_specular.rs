//! Fresnel-Modulated Specular Reflection and Transmission

#![allow(dead_code)]

use super::*;
use crate::material::*;
use std::fmt;

/// BRDF for physically plausible specular reflection and transmission.
#[derive(Clone)]
pub struct FresnelSpecular {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Spectrum used to scale the reflected colour.
    r: Spectrum,

    /// Spectrum used to scale the transmitted colour.
    t: Spectrum,

    /// Index of refraction above the surface (same side as surface normal).
    eta_a: Float,

    /// Index of refraction below the surface (opposite side as surface normal).
    eta_b: Float,

    /// Indicates whether incident ray started from a light source or from camera.
    mode: TransportMode,
}

impl FresnelSpecular {
    /// Creates a new instance of `FresnelSpecular`.
    ///
    /// * `fresnel` - Fresnel interface for dielectrics and conductors.
    /// * `r`       - Spectrum used to scale the reflected colour.
    /// * `t`       - Spectrum used to scale the transmitted colour.
    /// * `eta_a`   - Index of refraction above the surface (same side as surface
    ///               normal).
    /// * `eta_b`   - Index of refraction below the surface (opposite side as surface
    ///               normal).
    /// * `mode`    - Indicates whether incident ray started from a light source
    ///               or from camera.
    pub fn new(r: Spectrum, t: Spectrum, eta_a: Float, eta_b: Float, mode: TransportMode) -> BxDF {
        let model = Self {
            bxdf_type: BxDFType::BSDF_REFLECTION
                | BxDFType::BSDF_TRANSMISSION
                | BxDFType::BSDF_SPECULAR,
            r,
            t,
            eta_a,
            eta_b,
            mode,
        };
        BxDF::FresnelSpecular(model)
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
    /// * `wo`           - Outgoing direction.
    /// * `u`            - The 2D uniform random values.
    pub fn sample_f(&self, wo: &Vector3f, u: &Point2f) -> BxDFSample {
        let f = fr_dielectric(cos_theta(wo), self.eta_a, self.eta_b);

        if u[0] < f {
            // Compute specular reflection for `FresnelSpecular`.
            // Compute perfect specular reflection direction.
            let wi = Vector3f::new(-wo.x, -wo.y, wo.z);
            let sampled_type = BxDFType::BSDF_SPECULAR | BxDFType::BSDF_REFLECTION;
            let pdf = f;
            BxDFSample::new(f * self.r / abs_cos_theta(&wi), pdf, wi, sampled_type)
        } else {
            // Compute specular transmission for `FresnelSpecular`.
            // Figure out which `eta` is incident and which is transmitted.
            let entering = cos_theta(wo) > 0.0;
            let eta_i = if entering { self.eta_a } else { self.eta_b };
            let eta_t = if entering { self.eta_b } else { self.eta_a };

            // Compute ray direction for specular transmission.
            let sampled_type = BxDFType::BSDF_SPECULAR | BxDFType::BSDF_TRANSMISSION;
            if let Some(wi) = refract(
                wo,
                &Normal3f::new(0.0, 0.0, 1.0).face_forward(wo),
                eta_i / eta_t,
            ) {
                let mut ft = self.t * (1.0 - f);

                // Account for non-symmetry with transmission to different medium
                if self.mode == TransportMode::Radiance {
                    ft *= (eta_i * eta_i) / (eta_t * eta_t);
                }

                let pdf = 1.0 - f;
                BxDFSample::new(ft / abs_cos_theta(&wi), pdf, wi, sampled_type)
            } else {
                BxDFSample::from(sampled_type)
            }
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

impl fmt::Display for FresnelSpecular {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "FresnelSpecular {{ bxdf_type: {}, r: {}, t: {}, eta_a: {}, eta_b: {}, mode: {} }}",
            self.bxdf_type, self.r, self.t, self.eta_a, self.eta_b, self.mode,
        )
    }
}
