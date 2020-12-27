//! Scaled BxDF

#![allow(dead_code)]

use super::*;

/// BxDF scaling adapter scales a BxDF's contribution with a `Spectrum`.
#[derive(Clone)]
pub struct ScaledBxDF {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// The BxDF to scale.
    bxdf: ArcBxDF,

    /// Scaling value.
    scale: Spectrum,
}

impl BxDF for ScaledBxDF {
    /// Returns the BxDF type.
    fn get_type(&self) -> BxDFType {
        self.bxdf_type
    }

    /// Returns the value of the distribution function for the given pair of
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        self.scale * self.bxdf.f(wo, wi)
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    fn sample_f(&self, wo: &Vector3f, u: &Point2f) -> (Spectrum, Float, Vector3f, BxDFType) {
        let (sample, pdf, wi, bxdf_type) = self.bxdf.sample_f(wo, u);
        (self.scale * sample, pdf, wi, bxdf_type)
    }

    /// Evaluates the PDF for the sampling method.
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        self.bxdf.pdf(wo, wi)
    }

    /// Computes the hemispherical-directional reflectance function ρ.
    ///
    /// * `wo`      - Outgoing direction.
    /// * `samples` - Samples used b Monte Carlo algorithm.
    fn rho_hd(&self, wo: &Vector3f, samples: &[Point2f]) -> Spectrum {
        self.scale * self.bxdf.rho_hd(wo, samples)
    }

    /// Computes the hemispherical-hemispherical-directional reflectance function ρ.
    ///
    /// * `samples1` - Samples used b Monte Carlo algorithm.
    /// * `samples2` - Samples used b Monte Carlo algorithm.
    fn rho_hh(&self, samples1: &[Point2f], samples2: &[Point2f]) -> Spectrum {
        self.scale * self.bxdf.rho_hh(samples1, samples2)
    }
}
