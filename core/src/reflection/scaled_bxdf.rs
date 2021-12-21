//! Scaled BxDF

use super::*;
use bumpalo::Bump;

/// BxDF scaling adapter scales a BxDF's contribution with a `Spectrum`.
#[derive(Clone)]
pub struct ScaledBxDF {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// The BxDF to scale.
    bxdf: Box<BxDF>,

    /// Scaling value.
    scale: Spectrum,
}

impl ScaledBxDF {
    /// Create a new instance of `ScaledBxDF`.
    ///
    /// * `bxdf`  - The BxDF to scale.
    /// * `scale` - Scaling value.
    pub fn new(bxdf: BxDF, scale: Spectrum) -> Self {
        Self {
            bxdf_type: bxdf.get_type(),
            bxdf: Box::new(bxdf),
            scale,
        }
    }

    /// Allocate a new instance of `ScaledBxDF`.
    ///
    /// * `allocator` - The allocator.
    /// * `bxdf`      - The BxDF to scale.
    /// * `scale`     - Scaling value.
    pub fn alloc(allocator: &Bump, bxdf: BxDF, scale: Spectrum) -> BxDF {
        let model = allocator.alloc(Self::new(bxdf, scale)).to_owned();
        allocator.alloc(BxDF::ScaledBxDF(model)).to_owned()
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
