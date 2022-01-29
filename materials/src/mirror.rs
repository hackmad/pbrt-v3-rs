//! Mirror Material

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

/// Implements a simple mirror, modeled with perfect specular reflection.
pub struct MirrorMaterial {
    /// Reflectivity of the mirror.
    kr: ArcTexture<Spectrum>,

    /// Bump map.
    bump_map: Option<ArcTexture<Float>>,
}

impl MirrorMaterial {
    /// Create a new `MirrorMaterial`.
    ///
    /// * `kr`       - Reflectivity of the mirror.
    /// * `bump_map` - Optional bump map.
    pub fn new(kr: ArcTexture<Spectrum>, bump_map: Option<ArcTexture<Float>>) -> Self {
        Self { kr, bump_map }
    }
}

impl Material for MirrorMaterial {
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
            Material::bump(self, bump_map, si);
        }

        let bsdf = BSDF::alloc(arena, &si, None);

        // Evaluate textures for `MirrorMaterial` material and allocate BRDF.
        let r = self.kr.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        if !r.is_black() {
            let fresnel = FresnelNoOp::alloc(arena);
            let bxdf = SpecularReflection::alloc(arena, r, fresnel);
            bsdf.add(bxdf);
        }

        si.bsdf = Some(bsdf);
    }
}

impl From<&TextureParams> for MirrorMaterial {
    /// Create a mirror material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let kr = tp.get_spectrum_texture_or_else("Kr", Spectrum::new(0.9), |v| {
            Arc::new(ConstantTexture::new(v))
        });

        let bump_map = tp.get_float_texture_or_none("bumpmap");

        Self::new(kr, bump_map)
    }
}
