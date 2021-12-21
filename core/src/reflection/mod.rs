//! Reflection and surface scattering models

#![allow(dead_code)]
use crate::geometry::*;
use crate::pbrt::*;
use crate::sampling::*;
use crate::spectrum::*;
use std::sync::Arc;

mod bsdf;
mod bsdf_reader;
mod bxdf_sample;
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
pub use bsdf::*;
pub use bxdf_sample::*;
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

/// BxDF for BRDFs and BTDFs.
#[derive(Clone)]
pub enum BxDF {
    FourierBSDF(FourierBSDF),
    FresnelBlend(FresnelBlend),
    FresnelSpecular(FresnelSpecular),
    LambertianReflection(LambertianReflection),
    MicrofacetReflection(MicrofacetReflection),
    MicrofacetTransmission(MicrofacetTransmission),
    OrenNayar(OrenNayar),
    ScaledBxDF(ScaledBxDF),
    SpecularReflection(SpecularReflection),
    SpecularTransmission(SpecularTransmission),
}

impl BxDF {
    /// Returns the BxDF type.
    pub fn get_type(&self) -> BxDFType {
        match self {
            BxDF::FourierBSDF(bxdf) => bxdf.get_type(),
            BxDF::FresnelBlend(bxdf) => bxdf.get_type(),
            BxDF::FresnelSpecular(bxdf) => bxdf.get_type(),
            BxDF::LambertianReflection(bxdf) => bxdf.get_type(),
            BxDF::MicrofacetReflection(bxdf) => bxdf.get_type(),
            BxDF::MicrofacetTransmission(bxdf) => bxdf.get_type(),
            BxDF::OrenNayar(bxdf) => bxdf.get_type(),
            BxDF::ScaledBxDF(bxdf) => bxdf.get_type(),
            BxDF::SpecularReflection(bxdf) => bxdf.get_type(),
            BxDF::SpecularTransmission(bxdf) => bxdf.get_type(),
        }
    }

    /// Returns true if the reflection models match.
    ///
    /// * `t` - The reflection model to compare.
    pub fn matches_flags(&self, t: BxDFType) -> bool {
        let bxdf_type = self.get_type();
        bxdf_type & t == bxdf_type
    }

    /// Returns the value of the distribution function for the given pair of
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        match self {
            BxDF::FourierBSDF(bxdf) => bxdf.f(wo, wi),
            BxDF::FresnelBlend(bxdf) => bxdf.f(wo, wi),
            BxDF::FresnelSpecular(bxdf) => bxdf.f(wo, wi),
            BxDF::LambertianReflection(bxdf) => bxdf.f(wo, wi),
            BxDF::MicrofacetReflection(bxdf) => bxdf.f(wo, wi),
            BxDF::MicrofacetTransmission(bxdf) => bxdf.f(wo, wi),
            BxDF::OrenNayar(bxdf) => bxdf.f(wo, wi),
            BxDF::ScaledBxDF(bxdf) => bxdf.f(wo, wi),
            BxDF::SpecularReflection(bxdf) => bxdf.f(wo, wi),
            BxDF::SpecularTransmission(bxdf) => bxdf.f(wo, wi),
        }
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    pub fn sample_f(&self, wo: &Vector3f, u: &Point2f) -> BxDFSample {
        match self {
            BxDF::FourierBSDF(bxdf) => bxdf.sample_f(wo, u),
            BxDF::FresnelBlend(bxdf) => bxdf.sample_f(wo, u),
            BxDF::FresnelSpecular(bxdf) => bxdf.sample_f(wo, u),
            BxDF::MicrofacetReflection(bxdf) => bxdf.sample_f(wo, u),
            BxDF::MicrofacetTransmission(bxdf) => bxdf.sample_f(wo, u),
            BxDF::ScaledBxDF(bxdf) => bxdf.sample_f(wo, u),
            BxDF::SpecularReflection(bxdf) => bxdf.sample_f(wo, u),
            BxDF::SpecularTransmission(bxdf) => bxdf.sample_f(wo, u),
            _ => {
                // Cosine-sample the hemisphere, flipping the direction if necessary.
                let mut wi = cosine_sample_hemisphere(u);
                if wo.z < 0.0 {
                    wi.z *= -1.0;
                }
                let pdf = self.pdf(wo, &wi);
                BxDFSample::new(self.f(wo, &wi), pdf, wi, BxDFType::BSDF_NONE) // self.get_type() ???
            }
        }
    }

    /// Evaluates the PDF for the sampling method. Default is based on the
    /// cosine-weighted sampling in `BxDF::sample_f()` default implementation.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        match self {
            BxDF::FourierBSDF(bxdf) => bxdf.pdf(wo, wi),
            BxDF::FresnelBlend(bxdf) => bxdf.pdf(wo, wi),
            BxDF::FresnelSpecular(bxdf) => bxdf.pdf(wo, wi),
            BxDF::MicrofacetReflection(bxdf) => bxdf.pdf(wo, wi),
            BxDF::MicrofacetTransmission(bxdf) => bxdf.pdf(wo, wi),
            BxDF::ScaledBxDF(bxdf) => bxdf.pdf(wo, wi),
            BxDF::SpecularReflection(bxdf) => bxdf.pdf(wo, wi),
            BxDF::SpecularTransmission(bxdf) => bxdf.pdf(wo, wi),
            _ => {
                if same_hemisphere(wo, wi) {
                    abs_cos_theta(wi) * INV_PI
                } else {
                    0.0
                }
            }
        }
    }

    /// Computes the hemispherical-directional reflectance function ρ.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - Samples used by Monte Carlo algorithm.
    pub fn rho_hd(&self, wo: &Vector3f, u: &[Point2f]) -> Spectrum {
        match self {
            BxDF::LambertianReflection(bxdf) => bxdf.rho_hd(wo, u),
            BxDF::ScaledBxDF(bxdf) => bxdf.rho_hd(wo, u),
            _ => {
                let mut r = Spectrum::new(0.0);
                for s in u {
                    // Estimate one term of `rho_hd`.
                    let sample = self.sample_f(wo, s);
                    if sample.pdf > 0.0 {
                        r += sample.f * abs_cos_theta(&sample.wi) / sample.pdf;
                    }
                }
                r / u.len() as Float
            }
        }
    }

    /// Computes the hemispherical-hemispherical-directional reflectance function ρ.
    ///
    /// * `u1` - Samples used b Monte Carlo algorithm.
    /// * `u2` - Samples used b Monte Carlo algorithm.
    pub fn rho_hh(&self, u1: &[Point2f], u2: &[Point2f]) -> Spectrum {
        match self {
            BxDF::LambertianReflection(bxdf) => bxdf.rho_hh(u1, u2),
            BxDF::ScaledBxDF(bxdf) => bxdf.rho_hh(u1, u2),
            _ => {
                assert!(u1.len() == u2.len());

                let mut r = Spectrum::new(0.0);
                for (s1, s2) in u1.iter().zip(u2.iter()) {
                    // Estimate one term of `rho_hh`.
                    let wo = uniform_sample_hemisphere(s1);
                    let pdfo = uniform_hemisphere_pdf();
                    let sample = self.sample_f(&wo, s2);
                    let pdfi = sample.pdf;
                    if pdfi > 0.0 {
                        r += sample.f * abs_cos_theta(&sample.wi) * abs_cos_theta(&wo)
                            / (pdfo * pdfi);
                    }
                }
                r / (PI * u1.len() as Float)
            }
        }
    }
}

macro_rules! bxdf_from {
    ($struct: ty, $enum: ident) => {
        impl From<$struct> for BxDF {
            /// Wraps $struct in BxDF::$enum.
            fn from(bxdf: $struct) -> Self {
                Self::$enum(bxdf)
            }
        }
    };
}
bxdf_from!(FourierBSDF, FourierBSDF);
bxdf_from!(FresnelBlend, FresnelBlend);
bxdf_from!(FresnelSpecular, FresnelSpecular);
bxdf_from!(LambertianReflection, LambertianReflection);
bxdf_from!(MicrofacetReflection, MicrofacetReflection);
bxdf_from!(MicrofacetTransmission, MicrofacetTransmission);
bxdf_from!(OrenNayar, OrenNayar);
bxdf_from!(ScaledBxDF, ScaledBxDF);
bxdf_from!(SpecularReflection, SpecularReflection);
bxdf_from!(SpecularTransmission, SpecularTransmission);
