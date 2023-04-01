//! Medium Interactions

#![allow(dead_code)]
use super::Hit;
use crate::geometry::*;
use crate::medium::*;
use crate::pbrt::*;

/// MediumInteraction represents an interaction point in a scattering medium.
#[derive(Clone)]
pub struct MediumInteraction {
    /// The common interaction data.
    pub hit: Hit,

    /// The phase function.
    pub phase: PhaseFunction,
}

impl MediumInteraction {
    /// Create a new medium interaction.
    ///
    /// * `p`      - The point of interaction.
    /// * `wo`     - The negative ray direction (outgoing direction used
    ///              hen computing lighting at points).
    /// * `time`   - Time when interaction occurred.
    /// * `medium` - The medium.
    /// * `phase`  - The phase function.
    pub fn new(p: Point3f, wo: Vector3f, time: Float, medium: Option<ArcMedium>, phase: PhaseFunction) -> Self {
        Self {
            hit: Hit::new(
                p,
                time,
                Vector3f::ZERO,
                wo,
                Normal3f::ZERO,
                medium.map(MediumInterface::from),
            ),
            phase,
        }
    }

    /// Spawn's a new ray in the given direction.
    ///
    /// * `d` - The new direction.
    pub fn spawn_ray(&self, d: &Vector3f) -> Ray {
        self.hit.spawn_ray(d)
    }

    /// Spawn's a new ray towards another point.
    ///
    /// * `p` - The target point.
    pub fn spawn_ray_to_point(&self, p: &Point3f) -> Ray {
        self.hit.spawn_ray_to_point(p)
    }

    /// Spawn's a new ray towards another interaction.
    ///
    /// * `hit` - The interaction.
    pub fn spawn_ray_to_hit(&self, hit: &Hit) -> Ray {
        self.hit.spawn_ray_to_hit(hit)
    }
}
