//! Translucent Material

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

/// Implements a material that describes diffuse and glossy specular reflection
/// and transmission through the surface.
pub struct TranslucentMaterial {
    /// Coefficient of diffuse reflection and transmission.
    kd: ArcTexture<Spectrum>,

    /// Coefficient of glossy reflection and transmission.
    ks: ArcTexture<Spectrum>,

    /// Roughness.
    roughness: ArcTexture<Float>,

    /// Fraction of reflected light.
    reflect: ArcTexture<Spectrum>,

    /// Fraction of transmitted light.
    transmit: ArcTexture<Spectrum>,

    /// Bump map.
    bump_map: Option<ArcTexture<Float>>,

    /// Remap roughness value to [0, 1] where higher values represent larger
    /// highlights. If this is `false`, use the microfacet distributions `alpha`
    /// parameter.
    remap_roughness: bool,
}

impl TranslucentMaterial {
    /// Create a new `TranslucentMaterial`.
    ///
    /// * `kd`              - Coefficient of diffuse reflection and transmission.
    /// * `ks`              - Coefficient of glossy reflection and transmission.
    /// * `roughness`       - Roughness.
    /// * `reflect`         - Fraction of reflected light.
    /// * `transmit`        - Fraction of transmitted light.
    /// * `remap_roughness` - Remap roughness value to [0, 1] where higher values
    ///                       represent larger highlights. If this is `false`,
    ///                       use the microfacet distributions `alpha` parameter.
    /// * `bump_map`        - Optional bump map.
    pub fn new(
        kd: ArcTexture<Spectrum>,
        ks: ArcTexture<Spectrum>,
        roughness: ArcTexture<Float>,
        reflect: ArcTexture<Spectrum>,
        transmit: ArcTexture<Spectrum>,
        remap_roughness: bool,
        bump_map: Option<ArcTexture<Float>>,
    ) -> Self {
        Self {
            kd,
            ks,
            roughness,
            reflect,
            transmit,
            remap_roughness,
            bump_map,
        }
    }
}

impl Material for TranslucentMaterial {
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
        mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = &self.bump_map {
            Material::bump(self, bump_map, si);
        }

        const ETA: Float = 1.5;
        let bsdf = BSDF::alloc(arena, &si.hit, &si.shading, Some(ETA));

        // Evaluate textures for `TranslucentMaterial` material and allocate BRDF.
        let r = self
            .reflect
            .evaluate(&si.hit, &si.uv, &si.der)
            .clamp_default();
        let t = self
            .transmit
            .evaluate(&si.hit, &si.uv, &si.der)
            .clamp_default();
        if r.is_black() && t.is_black() {
            return;
        }

        let kd = self.kd.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        if !kd.is_black() {
            if !r.is_black() {
                bsdf.add(LambertianReflection::alloc(arena, r * kd));
            }
            if !t.is_black() {
                bsdf.add(LambertianTransmission::alloc(arena, t * kd));
            }
        }

        let ks = self.ks.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        if !ks.is_black() && (!r.is_black() || !t.is_black()) {
            let mut rough = self.roughness.evaluate(&si.hit, &si.uv, &si.der);
            if self.remap_roughness {
                rough = TrowbridgeReitzDistribution::roughness_to_alpha(rough);
            }

            if !r.is_black() {
                let distrib = TrowbridgeReitzDistribution::alloc(arena, rough, rough, true);
                let fresnel = FresnelDielectric::alloc(arena, 1.0, ETA);
                bsdf.add(MicrofacetReflection::alloc(arena, r * ks, distrib, fresnel));
            }

            if !t.is_black() {
                let distrib = TrowbridgeReitzDistribution::alloc(arena, rough, rough, true);
                let fresnel = FresnelDielectric::alloc(arena, 1.0, ETA);
                bsdf.add(MicrofacetTransmission::alloc(
                    arena,
                    t * ks,
                    distrib,
                    fresnel,
                    1.0,
                    ETA,
                    mode,
                ));
            }
        }

        si.bsdf = Some(bsdf);
    }
}

impl From<&TextureParams> for TranslucentMaterial {
    /// Create a translucent material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let kd = tp.get_spectrum_texture_or_else("Kd", Spectrum::new(0.25), |v| {
            Arc::new(ConstantTexture::new(v))
        });
        let ks = tp.get_spectrum_texture_or_else("Ks", Spectrum::new(0.25), |v| {
            Arc::new(ConstantTexture::new(v))
        });
        let reflect = tp.get_spectrum_texture_or_else("reflect", Spectrum::new(0.5), |v| {
            Arc::new(ConstantTexture::new(v))
        });
        let transmit = tp.get_spectrum_texture_or_else("transmit", Spectrum::new(0.5), |v| {
            Arc::new(ConstantTexture::new(v))
        });
        let roughness =
            tp.get_float_texture_or_else("roughness", 0.1, |v| Arc::new(ConstantTexture::new(v)));

        let bump_map = tp.get_float_texture_or_none("bumpmap");
        let remap_roughness = tp.find_bool("remaproughness", true);

        Self::new(
            kd,
            ks,
            roughness,
            reflect,
            transmit,
            remap_roughness,
            bump_map,
        )
    }
}
