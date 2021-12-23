//! Plastic Material

use bumpalo::Bump;
use core::interaction::*;
use core::material::*;
use core::microfacet::*;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::*;
use core::spectrum::*;
use core::texture::*;
use std::rc::Rc;
use std::sync::Arc;
use textures::*;

/// Implements plastic material.
pub struct PlasticMaterial {
    /// Spectral diffuse reflection.
    kd: ArcTexture<Spectrum>,

    /// Spectral specular reflection.
    ks: ArcTexture<Spectrum>,

    /// Roughness.
    roughness: ArcTexture<Float>,

    /// Bump map.
    bump_map: Option<ArcTexture<Float>>,

    /// Remap roughness value to [0, 1] where higher values represent larger
    /// highlights. If this is `false`, use the microfacet distributions `alpha`
    /// parameter.
    remap_roughness: bool,
}

impl PlasticMaterial {
    /// Create a new `PlasticMaterial`.
    ///
    ///
    /// * `kd`              - Spectral diffuse reflection.
    /// * `ks`              - Spectral specular reflection.
    /// * `roughness`       - Roughness.
    /// * `remap_roughness` - Remap roughness value to [0, 1] where higher values
    ///                       represent larger highlights. If this is `false`,
    ///                       use the microfacet distributions `alpha` parameter.
    /// * `bump_map`        - Optional bump map.
    pub fn new(
        kd: ArcTexture<Spectrum>,
        ks: ArcTexture<Spectrum>,
        roughness: ArcTexture<Float>,
        remap_roughness: bool,
        bump_map: Option<ArcTexture<Float>>,
    ) -> Self {
        Self {
            kd: Arc::clone(&kd),
            ks: Arc::clone(&ks),
            roughness: Arc::clone(&roughness),
            remap_roughness,
            bump_map: bump_map.clone(),
        }
    }
}

impl Material for PlasticMaterial {
    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// * `arena`                - The memory arena for allocations.
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode (ignored).
    /// * `allow_multiple_lobes` - Indicates whether the material should use
    ///                            BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs
    ///                            are available (ignored).
    fn compute_scattering_functions(
        &self,
        arena: &Bump,
        si: &mut SurfaceInteraction,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = &self.bump_map {
            Material::bump(self, Arc::clone(bump_map), si);
        }

        let mut bsdf = BSDF::alloc(arena, &si, None);

        // Initialize diffuse component of plastic material.
        let kd = self.kd.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        if !kd.is_black() {
            bsdf.add(LambertianReflection::alloc(arena, kd));
        }

        // Initialize specular component of plastic material.
        let ks = self.ks.evaluate(&si.hit, &si.uv, &si.der).clamp_default();
        if !ks.is_black() {
            let fresnel = FresnelDielectric::alloc(arena, 1.5, 1.0);

            // Create microfacet distribution for plastic material.
            let mut rough = self.roughness.evaluate(&si.hit, &si.uv, &si.der);
            if self.remap_roughness {
                rough = TrowbridgeReitzDistribution::roughness_to_alpha(rough);
            }
            let distrib = Rc::new(TrowbridgeReitzDistribution::alloc(
                arena, rough, rough, true,
            ));
            bsdf.add(MicrofacetReflection::alloc(arena, ks, distrib, fresnel));
        }

        si.bsdf = Some(bsdf);
    }
}

impl From<&TextureParams> for PlasticMaterial {
    /// Create a plastic material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let kd = tp.get_spectrum_texture_or_else("Kd", Spectrum::new(0.25), |v| {
            Arc::new(ConstantTexture::new(v))
        });

        let ks = tp.get_spectrum_texture_or_else("Ks", Spectrum::new(0.25), |v| {
            Arc::new(ConstantTexture::new(v))
        });

        let roughness =
            tp.get_float_texture_or_else("roughness", 0.1, |v| Arc::new(ConstantTexture::new(v)));

        let bump_map = tp.get_float_texture("bumpmap");
        let remap_roughness = tp.find_bool("remaproughness", true);

        Self::new(kd, ks, roughness, remap_roughness, bump_map)
    }
}
