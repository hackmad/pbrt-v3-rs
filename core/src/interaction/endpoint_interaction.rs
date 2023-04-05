//! Endpiont Interactions

use crate::camera::*;
use crate::geometry::*;
use crate::interaction::*;
use crate::light::*;
use std::sync::Arc;

/// Represents an interaction point used only by BDPT integrator.
#[derive(Clone)]
pub enum EndpointInteraction {
    /// Records the position of a path endpoint on the lens of the camera.
    Camera {
        /// The interaction point.
        hit: Hit,
        /// The camera.
        camera: ArcCamera,
    },

    /// Records the position of a path endpoint on a light source.
    Light {
        /// The interaction point.
        hit: Hit,
        /// The light source.
        light: Option<ArcLight>,
    },
}

impl EndpointInteraction {
    /// Create a camera endpoint interaction.
    ///
    /// * `hit`    - The hit point on camera lens.
    /// * `camera` - The camera.
    pub fn camera_from_hit(hit: Hit, camera: ArcCamera) -> Self {
        Self::Camera { hit, camera }
    }

    /// Create a camera endpoint interaction from a ray.
    ///
    /// * `ray`    - The ray starting on camera lens.
    /// * `camera` - The camera.
    pub fn camera_from_ray(ray: &Ray, camera: ArcCamera) -> Self {
        Self::Camera {
            hit: Hit::new(
                ray.o,
                ray.time,
                Vector3f::ZERO,
                Vector3f::ZERO,
                Normal3f::ZERO,
                ray.medium.as_ref().map(Arc::clone).map(MediumInterface::from),
            ),
            camera,
        }
    }

    /// Create a light endpoint interaction.
    ///
    /// * `hit`    - The hit point on a light source.
    /// * `light`  - The light source.
    pub fn light_from_hit(hit: Hit, light: Option<ArcLight>) -> Self {
        Self::Light { hit, light }
    }

    /// Create a light endpoint interaction from a ray.
    ///
    /// * `ray` - The ray starting point on a light source.
    pub fn light_from_ray(ray: &Ray) -> Self {
        let n = Normal3f::from(-ray.d);
        let hit = Hit::new(
            ray.at(1.0),
            ray.time,
            Vector3f::ZERO,
            Vector3f::ZERO,
            n,
            ray.medium.as_ref().map(Arc::clone).map(MediumInterface::from),
        );
        Self::Light { hit, light: None }
    }

    /// Create a light endpoint interaction from a ray and normal.
    ///
    /// * `ray`    - The ray starting point on a light source.
    /// * `normal` - The normal.
    /// * `light`  - The light source.
    pub fn light_from_ray_and_normal(ray: &Ray, nl: Normal3f, light: ArcLight) -> Self {
        Self::Light {
            hit: Hit::new(
                ray.o,
                ray.time,
                Vector3f::ZERO,
                Vector3f::ZERO,
                nl,
                ray.medium.as_ref().map(Arc::clone).map(MediumInterface::from),
            ),
            light: Some(light),
        }
    }

    /// Returns the hit point.
    pub fn hit(&self) -> &Hit {
        match self {
            Self::Camera { hit, .. } => &hit,
            Self::Light { hit, .. } => &hit,
        }
    }

    /// Spawn's a new ray in the given direction.
    ///
    /// * `d` - The new direction.
    pub fn spawn_ray(&self, d: &Vector3f) -> Ray {
        match self {
            Self::Camera { hit, .. } => hit.spawn_ray(d),
            Self::Light { hit, .. } => hit.spawn_ray(d),
        }
    }

    /// Spawn's a new ray towards another point.
    ///
    /// * `p` - The target point.
    pub fn spawn_ray_to_point(&self, p: &Point3f) -> Ray {
        match self {
            Self::Camera { hit, .. } => hit.spawn_ray_to_point(p),
            Self::Light { hit, .. } => hit.spawn_ray_to_point(p),
        }
    }

    /// Spawn's a new ray towards another interaction.
    ///
    /// * `h` - The interaction.
    pub fn spawn_ray_to_hit(&self, h: &Hit) -> Ray {
        match self {
            Self::Camera { hit, .. } => hit.spawn_ray_to_hit(h),
            Self::Light { hit, .. } => hit.spawn_ray_to_hit(h),
        }
    }
}
