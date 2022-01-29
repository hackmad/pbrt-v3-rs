//! Lambertian Reflection

#![allow(dead_code)]

use super::*;
use bumpalo::Bump;

/// BRDF for the Lambertian model for perfect transmissive surfaces that scatters
/// incident illumination equally through a surface in all directions.
pub struct LambertianTransmission {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Transmission spectrum which gives the fraction of incident light that
    /// is scattered through the surface.
    t: Spectrum,
}

impl LambertianTransmission {
    /// Allocate a new instance of `LambertianTransmission`.
    ///
    /// * `arena` - The arena for memory allocations.
    /// * `t`     - Transmission spectrum which gives the fraction of incident
    ///             light that is scattered through the surface.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc<'arena>(arena: &'arena Bump, t: Spectrum) -> &'arena mut BxDF {
        let model = arena.alloc(Self {
            bxdf_type: BxDFType::BSDF_TRANSMISSION | BxDFType::BSDF_DIFFUSE,
            t,
        });
        arena.alloc(BxDF::LambertianTransmission(model))
    }

    /// Clone into a newly allocated a new instance of `LambertianTransmission`.
    ///
    /// * `arena` - The arena for memory allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc<'arena>(&self, arena: &'arena Bump) -> &'arena mut BxDF<'arena> {
        let model = arena.alloc(Self {
            bxdf_type: self.bxdf_type,
            t: self.t,
        });
        arena.alloc(BxDF::LambertianTransmission(model))
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
        self.t * INV_PI
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    pub fn sample_f(&self, wo: &Vector3f, u: &Point2f) -> BxDFSample {
        // Cosine-sample the hemisphere, flipping the direction if necessary.
        let mut wi = cosine_sample_hemisphere(u);
        if wo.z > 0.0 {
            wi.z *= -1.0;
        }
        let pdf = self.pdf(wo, &wi);
        BxDFSample::new(self.f(wo, &wi), pdf, wi, self.get_type())
    }

    /// Evaluates the PDF for the sampling method. Default is based on the
    /// cosine-weighted sampling in `BxDF::sample_f()` default implementation.
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if !same_hemisphere(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0
        }
    }

    /// Computes the hemispherical-directional reflectance function ρ.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - Samples used by Monte Carlo algorithm.
    pub fn rho_hd(&self, _wo: &Vector3f, _u: &[Point2f]) -> Spectrum {
        self.t
    }

    /// Computes the hemispherical-hemispherical-directional reflectance function ρ.
    ///
    /// * `u1` - Samples used b Monte Carlo algorithm.
    /// * `u2` - Samples used b Monte Carlo algorithm.
    pub fn rho_hh(&self, u1: &[Point2f], u2: &[Point2f]) -> Spectrum {
        assert!(u1.len() == u2.len());
        self.t
    }
}
