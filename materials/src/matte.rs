//! Matte Material

use bumpalo::Bump;
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
    /// * `arena`                - The arena for memory allocations.
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode (ignored).
    /// * `allow_multiple_lobes` - Indicates whether the material should use
    ///                            BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs
    ///                            are available (ignored).
    fn compute_scattering_functions<'primtive, 'arena>(
        &self,
        arena: &'arena Bump,
        si: &mut SurfaceInteraction<'primtive, 'arena>,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = &self.bump_map {
            Material::bump(self, Arc::clone(bump_map), si);
        }

        let bsdf = BSDF::new(arena, &si, None);

        // Evaluate textures for `MatteMaterial` material and allocate BRDF
        let r = self.kd.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        let sig = clamp(self.sigma.evaluate(&si.hit, &si.uv, &si.der), 0.0, 90.0);
        if !r.is_black() {
            let bxdf = if sig == 0.0 {
                LambertianReflection::new(arena, r)
            } else {
                OrenNayar::new(arena, r, sig)
            };
            bsdf.add(bxdf);
        }

        si.bsdf = Some(bsdf);
    }
}

impl From<&TextureParams> for MatteMaterial {
    /// Create a matte material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let kd = tp.get_spectrum_texture_or_else("Kd", Spectrum::new(0.5), |v| {
            Arc::new(ConstantTexture::new(v))
        });

        let sigma =
            tp.get_float_texture_or_else("sigma", 0.0, |v| Arc::new(ConstantTexture::new(v)));

        let bump_map = tp.get_float_texture("bumpmap");

        Self::new(kd, sigma, bump_map)
    }
}
