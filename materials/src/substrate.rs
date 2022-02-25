//! Substrate Material

use bumpalo::Bump;
use core::interaction::*;
use core::material::*;
use core::microfacet::TrowbridgeReitzDistribution;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::*;
use core::spectrum::*;
use core::texture::*;
use std::sync::Arc;
use textures::*;

/// Implements a layered model that varies between glossy specular and diffuse
/// reflection depending on the viewing angle (based on the FresnelBlend BRDF).
pub struct SubstrateMaterial {
    /// Coefficient of diffuse reflection.
    kd: ArcTexture<Spectrum>,

    /// Coefficient of specular reflection.
    ks: ArcTexture<Spectrum>,

    /// Microfacet roughness in the u direction. If zero, perfect specular
    /// reflection is modeled.
    u_roughness: ArcTexture<Float>,

    /// Microfacet roughness in the v direction. If zero, perfect specular
    /// reflection is modeled.
    v_roughness: ArcTexture<Float>,

    /// If true, roughness values are expected to be in the range [0,1], and are
    /// remapped to microfacet distribution function parameter values that range
    /// from near-perfect-specular at 0 to very rough at 1. Otherwise the
    /// roughness parameters are used directly for the alpha parameters of the
    /// microfacet distribution function.
    remap_roughness: bool,

    /// Bump map.
    bump_map: Option<ArcTexture<Float>>,
}

impl SubstrateMaterial {
    /// Create a new `SubstrateMaterial`.
    ///
    /// * `kd`              - Coefficient of diffuse reflection.
    /// * `ks`              - Coefficient of specular reflection.
    /// * `u_roughness`     - Microfacet roughness in the u direction. If zero,
    ///                       perfect specular reflection is modeled.
    /// * `v_roughness`     - Microfacet roughness in the v direction. If zero,
    ///                       perfect specular reflection is modeled.
    /// * `remap_roughness` - If true, roughness values are expected to be in
    ///                       the range [0,1], and are remapped to microfacet
    ///                       distribution function parameter values that range
    ///                       from near-perfect-specular at 0 to very rough at 1.
    ///                       Otherwise the roughness parameters are used directly
    ///                       for the alpha parameters of the microfacet distribution
    ///                       function.
    /// * `bump_map`        - Optional bump map.
    pub fn new(
        kd: ArcTexture<Spectrum>,
        ks: ArcTexture<Spectrum>,
        u_roughness: ArcTexture<Float>,
        v_roughness: ArcTexture<Float>,
        remap_roughness: bool,
        bump_map: Option<ArcTexture<Float>>,
    ) -> Self {
        Self {
            kd,
            ks,
            u_roughness,
            v_roughness,
            remap_roughness,
            bump_map,
        }
    }
}

impl Material for SubstrateMaterial {
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
    fn compute_scattering_functions<'scene, 'arena>(
        &self,
        arena: &'arena Bump,
        si: &mut SurfaceInteraction<'scene, 'arena>,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = &self.bump_map {
            Material::bump(self, bump_map, si);
        }

        let bsdf = BSDF::alloc(arena, &si.hit, &si.shading, None);

        // Evaluate textures for `SubstrateMaterial` material and allocate BRDF.
        let d = self.kd.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        let s = self.ks.evaluate(&si.hit, &si.uv, &si.der).clamp_default();

        let mut urough = self.u_roughness.evaluate(&si.hit, &si.uv, &si.der);
        let mut vrough = self.v_roughness.evaluate(&si.hit, &si.uv, &si.der);

        if !d.is_black() || !s.is_black() {
            if self.remap_roughness {
                urough = TrowbridgeReitzDistribution::roughness_to_alpha(urough);
                vrough = TrowbridgeReitzDistribution::roughness_to_alpha(vrough);
            }
            let distrib = TrowbridgeReitzDistribution::alloc(arena, urough, vrough, true);
            bsdf.add(FresnelBlend::alloc(arena, d, s, distrib));
        }

        si.bsdf = Some(bsdf);
    }
}

impl From<&TextureParams> for SubstrateMaterial {
    /// Create a substrate material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let kd = tp.get_spectrum_texture_or_else("Kd", Spectrum::new(0.5), |v| {
            Arc::new(ConstantTexture::new(v))
        });

        let ks = tp.get_spectrum_texture_or_else("Ks", Spectrum::new(0.5), |v| {
            Arc::new(ConstantTexture::new(v))
        });

        let u_roughness =
            tp.get_float_texture_or_else("uroughness", 0.1, |v| Arc::new(ConstantTexture::new(v)));

        let v_roughness =
            tp.get_float_texture_or_else("vroughness", 0.1, |v| Arc::new(ConstantTexture::new(v)));

        let bump_map = tp.get_float_texture_or_none("bumpmap");
        let remap_roughness = tp.find_bool("remaproughness", true);

        Self::new(kd, ks, u_roughness, v_roughness, remap_roughness, bump_map)
    }
}
