//! Matte Material

use core::geometry::*;
use core::material::*;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::*;
use core::spectrum::*;
use core::texture::*;
use std::sync::Arc;
use textures::*;

/// Implements purely diffuse surfaces.
pub struct MatteMaterial {
    /// Spectral diffuse reflection.
    kd: ArcTexture<Spectrum>,

    /// Roughness.
    sigma: ArcTexture<Float>,

    /// Bump map.
    bump_map: Option<ArcTexture<Float>>,
}

impl MatteMaterial {
    /// Create a new `MatteMaterial`.
    ///
    ///
    /// * `kd`       - Spectral diffuse reflection.
    /// * `sigma`    - Roughness.
    /// * `bump_map` - Optional bump map.
    pub fn new(
        kd: ArcTexture<Spectrum>,
        sigma: ArcTexture<Float>,
        bump_map: Option<ArcTexture<Float>>,
    ) -> Self {
        Self {
            kd: Arc::clone(&kd),
            sigma: Arc::clone(&sigma),
            bump_map: bump_map.clone(),
        }
    }
}

impl Material for MatteMaterial {
    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode (ignored).
    /// * `allow_multiple_lobes` - Indicates whether the material should use
    ///                            BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs
    ///                            are available (ignored).
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = self.bump_map.clone() {
            Material::bump(self, bump_map, si);
        }

        let mut bsdf = BSDF::new(&si, None);

        // Evaluate textures for `MatteMaterial` material and allocate BRDF
        let r = self.kd.evaluate(si).clamp_default();
        let sig = clamp(self.sigma.evaluate(si), 0.0, 90.0);
        if !r.is_black() {
            if sig == 0.0 {
                bsdf.add(Arc::new(LambertianReflection::new(r)));
            } else {
                bsdf.add(Arc::new(OrenNayar::new(r, sig)));
            }
        }

        si.bsdf = Some(bsdf);
    }
}

impl From<&TextureParams> for MatteMaterial {
    /// Create a matte material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let kd = tp
            .get_spectrum_texture_or_else("Kd", Arc::new(ConstantTexture::new(Spectrum::new(0.5))));
        let sigma = tp.get_float_texture_or_else("sigma", Arc::new(ConstantTexture::new(0.0)));
        let bump_map = tp.get_float_texture("bumpmap");
        Self::new(kd, sigma, bump_map)
    }
}
