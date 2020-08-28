//! Material

#![allow(dead_code)]
use super::geometry::SurfaceInteraction;
use super::pbrt::Float;
use super::texture::ArcTexture;
use std::sync::Arc;

// TransportMode enumeration.
pub enum TransportMode {
    Radiance,
    Importance,
}

/// Material trait provides common behavior.
pub trait Material {
    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode.
    /// * `allow_multiple_lobes` - Indicates whether the material should use
    ///                            BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs
    ///                            are available.
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    );

    /// Update the normal at the surface interaction using a bump map.
    ///
    /// * `d`  - Bump map.
    /// * `si` - Surface interaction.
    fn bump(&self, d: ArcTexture<Float>, si: &mut SurfaceInteraction);
}

/// Atomic reference counted `Material`.
pub type ArcMaterial = Arc<dyn Material + Send + Sync>;
