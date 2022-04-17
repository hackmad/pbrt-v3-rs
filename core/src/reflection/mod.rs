//! Reflection and surface scattering models

#![allow(dead_code)]
use crate::geometry::*;
use crate::pbrt::*;
use crate::sampling::*;
use crate::spectrum::*;
use bumpalo::Bump;
use std::fmt;

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
mod lambertian_transmission;
mod microfacet_reflection;
mod microfacet_transmission;
mod oren_nayar;
mod scaled_bxdf;
mod specular_reflection;
mod specular_transmission;
mod tabulated_bssrdf;

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
pub use lambertian_transmission::*;
pub use microfacet_reflection::*;
pub use microfacet_transmission::*;
pub use oren_nayar::*;
pub use scaled_bxdf::*;
pub use specular_reflection::*;
pub use specular_transmission::*;
pub use tabulated_bssrdf::*;

/// BxDF for BRDFs and BTDFs.
///
/// NOTES:
///
/// The PBRT source code uses a `SeparableBSDFAdapter`. We bypass that by
/// enumerating the BSSRDFs directly to avoid dealing with trait objects or
/// nesting enumerations which will add more boiler plate code.
pub enum BxDF<'arena> {
    FourierBSDF(&'arena mut FourierBSDF),
    FresnelBlend(&'arena mut FresnelBlend<'arena>),
    FresnelSpecular(&'arena mut FresnelSpecular),
    LambertianReflection(&'arena mut LambertianReflection),
    LambertianTransmission(&'arena mut LambertianTransmission),
    MicrofacetReflection(&'arena mut MicrofacetReflection<'arena>),
    MicrofacetTransmission(&'arena mut MicrofacetTransmission<'arena>),
    OrenNayar(&'arena mut OrenNayar),
    ScaledBxDF(&'arena mut ScaledBxDF<'arena>),
    SpecularReflection(&'arena mut SpecularReflection<'arena>),
    SpecularTransmission(&'arena mut SpecularTransmission<'arena>),
    TabulatedBSSRDF(&'arena mut TabulatedBSSRDF<'arena>),
}

impl<'arena> BxDF<'arena> {
    /// Clone into a newly allocated instance.
    ///
    /// * `arena` - The memory arena used for allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc(&self, arena: &'arena Bump) -> &'arena mut BxDF<'arena> {
        match self {
            BxDF::FourierBSDF(bxdf) => bxdf.clone_alloc(arena),
            BxDF::FresnelBlend(bxdf) => bxdf.clone_alloc(arena),
            BxDF::FresnelSpecular(bxdf) => bxdf.clone_alloc(arena),
            BxDF::LambertianReflection(bxdf) => bxdf.clone_alloc(arena),
            BxDF::LambertianTransmission(bxdf) => bxdf.clone_alloc(arena),
            BxDF::MicrofacetReflection(bxdf) => bxdf.clone_alloc(arena),
            BxDF::MicrofacetTransmission(bxdf) => bxdf.clone_alloc(arena),
            BxDF::OrenNayar(bxdf) => bxdf.clone_alloc(arena),
            BxDF::ScaledBxDF(bxdf) => bxdf.clone_alloc(arena),
            BxDF::SpecularReflection(bxdf) => bxdf.clone_alloc(arena),
            BxDF::SpecularTransmission(bxdf) => bxdf.clone_alloc(arena),
            BxDF::TabulatedBSSRDF(bxdf) => bxdf.clone_alloc(arena),
        }
    }

    /// Returns the BxDF type.
    pub fn get_type(&self) -> BxDFType {
        match self {
            BxDF::FourierBSDF(bxdf) => bxdf.get_type(),
            BxDF::FresnelBlend(bxdf) => bxdf.get_type(),
            BxDF::FresnelSpecular(bxdf) => bxdf.get_type(),
            BxDF::LambertianReflection(bxdf) => bxdf.get_type(),
            BxDF::LambertianTransmission(bxdf) => bxdf.get_type(),
            BxDF::MicrofacetReflection(bxdf) => bxdf.get_type(),
            BxDF::MicrofacetTransmission(bxdf) => bxdf.get_type(),
            BxDF::OrenNayar(bxdf) => bxdf.get_type(),
            BxDF::ScaledBxDF(bxdf) => bxdf.get_type(),
            BxDF::SpecularReflection(bxdf) => bxdf.get_type(),
            BxDF::SpecularTransmission(bxdf) => bxdf.get_type(),
            BxDF::TabulatedBSSRDF(bxdf) => bxdf.get_type(),
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
            BxDF::LambertianTransmission(bxdf) => bxdf.f(wo, wi),
            BxDF::MicrofacetReflection(bxdf) => bxdf.f(wo, wi),
            BxDF::MicrofacetTransmission(bxdf) => bxdf.f(wo, wi),
            BxDF::OrenNayar(bxdf) => bxdf.f(wo, wi),
            BxDF::ScaledBxDF(bxdf) => bxdf.f(wo, wi),
            BxDF::SpecularReflection(bxdf) => bxdf.f(wo, wi),
            BxDF::SpecularTransmission(bxdf) => bxdf.f(wo, wi),
            BxDF::TabulatedBSSRDF(bxdf) => bxdf.f(wo, wi),
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
            BxDF::LambertianTransmission(bxdf) => bxdf.sample_f(wo, u),
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
                BxDFSample::new(self.f(wo, &wi), pdf, wi, self.get_type())
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
            BxDF::LambertianTransmission(bxdf) => bxdf.pdf(wo, wi),
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
            BxDF::LambertianTransmission(bxdf) => bxdf.rho_hd(wo, u),
            BxDF::ScaledBxDF(bxdf) => bxdf.rho_hd(wo, u),
            _ => {
                let mut r = Spectrum::ZERO;
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
            BxDF::LambertianTransmission(bxdf) => bxdf.rho_hh(u1, u2),
            BxDF::ScaledBxDF(bxdf) => bxdf.rho_hh(u1, u2),
            _ => {
                assert!(u1.len() == u2.len());

                let mut r = Spectrum::ZERO;
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

impl<'arena> fmt::Display for BxDF<'arena> {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BxDF::FourierBSDF(bxdf) => write!(f, "BxDF {{ {} }}", bxdf),
            BxDF::FresnelBlend(bxdf) => write!(f, "BxDF {{ {} }}", bxdf),
            BxDF::FresnelSpecular(bxdf) => write!(f, "BxDF {{ {} }}", bxdf),
            BxDF::LambertianReflection(bxdf) => write!(f, "BxDF {{ {} }}", bxdf),
            BxDF::LambertianTransmission(bxdf) => write!(f, "BxDF {{ {} }}", bxdf),
            BxDF::MicrofacetReflection(bxdf) => write!(f, "BxDF {{ {} }}", bxdf),
            BxDF::MicrofacetTransmission(bxdf) => write!(f, "BxDF {{ {} }}", bxdf),
            BxDF::OrenNayar(bxdf) => write!(f, "BxDF {{ {} }}", bxdf),
            BxDF::ScaledBxDF(bxdf) => write!(f, "BxDF {{ {} }}", bxdf),
            BxDF::SpecularReflection(bxdf) => write!(f, "BxDF {{ {} }}", bxdf),
            BxDF::SpecularTransmission(bxdf) => write!(f, "BxDF {{ {} }}", bxdf),
            BxDF::TabulatedBSSRDF(bxdf) => write!(f, "BxDF {{ {} }}", bxdf),
        }
    }
}
