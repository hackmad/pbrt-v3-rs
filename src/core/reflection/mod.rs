//! Reflection and surface scattering models

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::sampling::*;
use crate::core::spectrum::*;
use std::sync::Arc;

mod bsdf_reader;
mod bxdf_type;
mod common;
mod fourier_bsdf;
mod fourier_bsdf_table;
mod fresnel;
mod fresnel_blend;
mod fresnel_specular;
mod lambertian_reflection;
mod microfacet_reflection;
mod microfacet_transmission;
mod oren_nayar;
mod scaled_bxdf;
mod specular_reflection;
mod specular_transmission;

// Re-export
pub use bxdf_type::*;
pub use common::*;
pub use fourier_bsdf::*;
pub use fourier_bsdf_table::*;
pub use fresnel::*;
pub use fresnel_blend::*;
pub use fresnel_specular::*;
pub use lambertian_reflection::*;
pub use microfacet_reflection::*;
pub use microfacet_transmission::*;
pub use oren_nayar::*;
pub use scaled_bxdf::*;
pub use specular_reflection::*;
pub use specular_transmission::*;

/// BxDF interface.
pub trait BxDF {
    /// Returns the BxDF type.
    fn get_type(&self) -> BxDFType;

    /// Returns `true` if the reflection model matches the given BxDF type.
    ///
    /// * `bxdf_type` - Type to compare.
    fn matches(&self, bxdf_type: BxDFType) -> bool {
        self.get_type() == bxdf_type
    }

    /// Returns the value of the distribution function for the given pair of
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum;

    /// Returns the value of the BxDF given the outgpoing direction.
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    fn sample_f(&self, wo: &Vector3f, u: &Point2f) -> (Spectrum, Float, Vector3f, BxDFType) {
        // Cosine-sample the hemisphere, flipping the direction if necessary.
        let mut wi = cosine_sample_hemisphere(u);
        if wo.z < 0.0 {
            wi.z *= -1.0;
        }
        let pdf = self.pdf(wo, &wi);
        (self.f(wo, &wi), pdf, wi, self.get_type())
    }

    /// Evaluates the PDF for the sampling method. Default is based on the
    /// cosine-weighted sampling in `BxDF::sample_f()` default implementation.
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if same_hemisphere(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0
        }
    }

    /// Computes the hemispherical-directional reflectance function ρ.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - Samples used by Monte Carlo algorithm.
    fn rho_hd(&self, wo: &Vector3f, u: &[Point2f]) -> Spectrum {
        let mut r = Spectrum::new(0.0);
        for s in u {
            // Estimate one term of `rho_hd`.
            let (f, pdf, wi, _sampled_type) = self.sample_f(wo, s);
            if pdf > 0.0 {
                r += f * abs_cos_theta(&wi) / pdf;
            }
        }
        r / u.len() as Float
    }

    /// Computes the hemispherical-hemispherical-directional reflectance function ρ.
    ///
    /// * `u1` - Samples used b Monte Carlo algorithm.
    /// * `u2` - Samples used b Monte Carlo algorithm.
    fn rho_hh(&self, u1: &[Point2f], u2: &[Point2f]) -> Spectrum {
        assert!(u1.len() == u2.len());

        let mut r = Spectrum::new(0.0);
        for (s1, s2) in u1.iter().zip(u2.iter()) {
            // Estimate one term of `rho_hh`.
            let wo = uniform_sample_hemisphere(s1);
            let pdfo = uniform_hemisphere_pdf();
            let (f, pdfi, wi, _sampled_type) = self.sample_f(&wo, s2);
            if pdfi > 0.0 {
                r += f * abs_cos_theta(&wi) * abs_cos_theta(&wo) / (pdfo * pdfi);
            }
        }
        r / (PI * u1.len() as Float)
    }
}

/// Atomic reference counted `BxDF`.
pub type ArcBxDF = Arc<dyn BxDF + Send + Sync>;

/// BSDF model.
pub struct BSDF {}

/// Atomic reference counted `BSDF`.
pub type ArcBSDF = Arc<BSDF>;
