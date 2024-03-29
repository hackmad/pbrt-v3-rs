//! Matte Material

use core::bssrdf::*;
use core::interaction::*;
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
    /// * `kd`       - Spectral diffuse reflection.
    /// * `sigma`    - Roughness.
    /// * `bump_map` - Optional bump map.
    pub fn new(kd: ArcTexture<Spectrum>, sigma: ArcTexture<Float>, bump_map: Option<ArcTexture<Float>>) -> Self {
        Self { kd, sigma, bump_map }
    }
}

impl Material for MatteMaterial {
    /// Initializes representations of the light-scattering properties of the material at the intersection point on the
    /// surface.
    ///
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode (ignored).
    /// * `allow_multiple_lobes` - Indicates whether the material should use BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs are available (ignored).
    /// * `bsdf`                 - The computed BSDF.
    /// * `bssrdf`               - The computed BSSSRDF.
    fn compute_scattering_functions<'scene>(
        &self,
        si: &mut SurfaceInteraction<'scene>,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
        bsdf: &mut Option<BSDF>,
        bssrdf: &mut Option<BSSRDF>,
    ) {
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = &self.bump_map {
            Material::bump(self, bump_map, si);
        }

        let mut result = BSDF::new(&si.hit, &si.shading, None);

        // Evaluate textures for `MatteMaterial` material and allocate BRDF.
        let r = self.kd.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        let sig = clamp(self.sigma.evaluate(&si.hit, &si.uv, &si.der), 0.0, 90.0);
        if !r.is_black() {
            let bxdf = if sig == 0.0 {
                LambertianReflection::new(r)
            } else {
                OrenNayar::new(r, sig)
            };
            result.add(bxdf);
        }

        *bsdf = Some(result);
        *bssrdf = None;
    }
}

impl From<&TextureParams> for MatteMaterial {
    /// Create a matte material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let kd = tp.get_spectrum_texture_or_else("Kd", Spectrum::new(0.5), |v| Arc::new(ConstantTexture::new(v)));

        let sigma = tp.get_float_texture_or_else("sigma", 0.0, |v| Arc::new(ConstantTexture::new(v)));

        let bump_map = tp.get_float_texture_or_none("bumpmap");

        Self::new(kd, sigma, bump_map)
    }
}
