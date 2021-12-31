//! Scaled BxDF

use super::*;
use bumpalo::Bump;

/// BxDF scaling adapter scales a BxDF's contribution with a `Spectrum`.
pub struct ScaledBxDF<'arena> {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// The BxDF to scale.
    ///
    /// NOTE: This is `&'arena mut BxDF` because it is allocated from a memory
    /// arena (bumpalo::Bump).
    bxdf: &'arena mut BxDF<'arena>,

    /// Scaling value.
    scale: Spectrum,
}

impl<'arena> ScaledBxDF<'arena> {
    /// Allocate a new instance of `ScaledBxDF`.
    ///
    /// * `arena` - The arena for memory allocations.
    /// * `bxdf`  - The BxDF to scale.
    /// * `scale` - Scaling value.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc(
        arena: &'arena Bump,
        bxdf: &'arena mut BxDF<'arena>,
        scale: Spectrum,
    ) -> &'arena mut BxDF<'arena> {
        let model = arena.alloc(Self {
            bxdf_type: bxdf.get_type(),
            bxdf,
            scale,
        });
        arena.alloc(BxDF::ScaledBxDF(model))
    }

    /// Clone into a newly allocated a new instance of `ScaledBxDF`.
    ///
    /// * `arena` - The arena for memory allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc(&self, arena: &'arena Bump) -> &'arena mut BxDF<'arena> {
        let bxdf = self.bxdf.clone_alloc(arena);
        let model = arena.alloc(Self {
            bxdf_type: self.bxdf_type,
            bxdf,
            scale: self.scale,
        });
        arena.alloc(BxDF::ScaledBxDF(model))
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
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        self.scale * self.bxdf.f(wo, wi)
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    pub fn sample_f(&self, wo: &Vector3f, u: &Point2f) -> BxDFSample {
        let mut sample = self.bxdf.sample_f(wo, u);
        sample.f = self.scale * sample.f;
        sample
    }

    /// Evaluates the PDF for the sampling method.
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        self.bxdf.pdf(wo, wi)
    }

    /// Computes the hemispherical-directional reflectance function ρ.
    ///
    /// * `wo`      - Outgoing direction.
    /// * `samples` - Samples used b Monte Carlo algorithm.
    pub fn rho_hd(&self, wo: &Vector3f, samples: &[Point2f]) -> Spectrum {
        self.scale * self.bxdf.rho_hd(wo, samples)
    }

    /// Computes the hemispherical-hemispherical-directional reflectance function ρ.
    ///
    /// * `samples1` - Samples used b Monte Carlo algorithm.
    /// * `samples2` - Samples used b Monte Carlo algorithm.
    pub fn rho_hh(&self, samples1: &[Point2f], samples2: &[Point2f]) -> Spectrum {
        self.scale * self.bxdf.rho_hh(samples1, samples2)
    }
}
