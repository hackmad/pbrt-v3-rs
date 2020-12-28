//! Matte Material

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::material::*;
use crate::core::pbrt::*;
use crate::core::reflection::*;
use crate::core::spectrum::*;
use crate::core::texture::*;
use std::sync::Arc;

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
            kd: kd.clone(),
            sigma: sigma.clone(),
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

        let mut bsdf = BSDF::new(&si.clone(), None);

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

        si.bsdf = Some(Arc::new(bsdf));
    }
}
