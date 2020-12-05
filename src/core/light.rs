//! Light

#![allow(dead_code)]
use super::geometry::ArcInteraction;
use super::sampler::ArcSampler;
use super::scene::Scene;
use super::spectrum::Spectrum;
use std::sync::Arc;

/// Light trait provides common behavior.
pub trait Light {}

/// Atomic reference counted `Light`.
pub type ArcLight = Arc<dyn Light + Send + Sync>;

/// AreaLight trait provides common behavior for area lights.
pub trait AreaLight: Light {}

/// Atomic reference counted `AreaLight`.
pub type ArcAreaLight = Arc<dyn AreaLight + Send + Sync>;

/// VisibilityTester allows lights to return a radiance value under the
/// assumption that the reference point and light source are mutually
/// visible.
#[derive(Clone)]
pub struct VisibilityTester {
    /// One endpoint of shadow ray.
    pub p0: ArcInteraction,

    /// Second endpoint of shadow ray.
    pub p1: ArcInteraction,
}

impl VisibilityTester {
    /// Create a new `VisibilityTester` for given endpoints of a shadow ray.
    ///
    /// * `p0` - One endpoint of shadow ray.
    /// * `p1` - Second endpoint of shadow ray.
    pub fn new(p0: ArcInteraction, p1: ArcInteraction) -> Self {
        Self {
            p0: p0.clone(),
            p1: p1.clone(),
        }
    }

    /// Traces a shadow ray between `p0` and `p1` through the scene and returns
    /// true if the points are visible to each other.
    ///
    /// * `scene` - The scene.
    pub fn unoccluded(&self, _scene: &Scene) -> bool {
        false
    }

    /// Computes the beam transmittance, the fraction of radiance transmitted
    /// along the segment between the two points. It accounts for both attenuation
    /// in participating media as well as any surfaces that block the ray completely.
    ///
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    pub fn tr(&self, _scene: &Scene, _sampler: ArcSampler) -> Spectrum {
        Spectrum::default()
    }
}
