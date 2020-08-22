//! Disks

#![allow(dead_code)]
use super::{
    bounds3, clamp, intersection, normal3, point2, point3, shape_data, surface_interaction,
    vector3, ArcTransform, Bounds3f, Float, Intersection, Ray, Shape, ShapeData, TWO_PI,
};
use std::sync::Arc;

/// A disk centered on the z-axis.
#[derive(Clone)]
pub struct Disk {
    /// Common shape data.
    pub data: ShapeData,

    /// Height where disk is located.
    pub height: Float,

    /// Radius of disk.
    pub radius: Float,

    /// Inner radius of disk to truncate center.
    pub inner_radius: Float,

    /// Maximum angle Φ to truncate disk.
    pub phi_max: Float,
}

/// Create a new disk centered on the z-axis.
///
/// * `object_to_world`     - The object to world transfomation.
/// * `world_to_object`     - The world to object transfomation.
/// * `reverse_orientation` - Indicates whether their surface normal directions
///                           should be reversed from the default
/// * `height`              - Height where disk is located.
/// * `radius`              - Radius of disk.
/// * `inner_radius`        - Inner radius of disk to truncate center.
/// * `phi_max`             - Maximum angle Φ to truncate disk.
pub fn disk(
    object_to_world: ArcTransform,
    world_to_object: ArcTransform,
    reverse_orientation: bool,
    height: Float,
    radius: Float,
    inner_radius: Float,
    phi_max: Float,
) -> Disk {
    Disk {
        height,
        radius,
        inner_radius,
        phi_max: clamp(phi_max, 0.0, 360.0).to_radians(),
        data: shape_data(
            object_to_world.clone(),
            Some(world_to_object.clone()),
            reverse_orientation,
        ),
    }
}

impl Shape for Disk {
    /// Returns the underlying shape data.
    fn get_data(&self) -> ShapeData {
        self.data.clone()
    }

    /// Returns a bounding box in the shapes object space.
    fn object_bound(&self) -> Bounds3f {
        bounds3(
            point3(-self.radius, -self.radius, self.height),
            point3(self.radius, self.radius, self.height),
        )
    }

    /// Returns geometric details if a ray intersects the shape intersection.
    /// If there is no intersection, `None` is returned.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests (not supported).
    fn intersect(&self, r: &Ray, _test_alpha_texture: bool) -> Option<Intersection> {
        // Transform ray to object space
        //
        // We could just use transform_ray() but there is minor adjustment in
        // it that adjusts t_max which is not in transform_ray_with_error().
        let (ray, _o_err, _d_err) = self
            .data
            .world_to_object
            .clone()
            .unwrap()
            .transform_ray_with_error(r);

        // Compute plane intersection for disk.

        // Reject disk intersections for rays parallel to the disk's plane
        if ray.d.z == 0.0 {
            return None;
        }
        let t_shape_hit = (self.height - ray.o.z) / ray.d.z;
        if t_shape_hit <= 0.0 || t_shape_hit >= ray.t_max {
            return None;
        }

        // See if hit point is inside disk radii and phimax
        let mut p_hit = ray.at(t_shape_hit);
        let dist2 = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
        if dist2 > self.radius * self.radius || dist2 < self.inner_radius * self.inner_radius {
            return None;
        }

        // Test disk phi value against phimax
        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += TWO_PI;
        }
        if phi > self.phi_max {
            return None;
        }

        // Find parametric representation of disk hit
        let u = phi / self.phi_max;
        let r_hit = dist2.sqrt();
        let v = (self.radius - r_hit) / (self.radius - self.inner_radius);
        let dpdu = vector3(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, 0.0);
        let dpdv = vector3(p_hit.x, p_hit.y, 0.0) * (self.inner_radius - self.radius) / r_hit;
        let dndu = normal3(0.0, 0.0, 0.0);
        let dndv = normal3(0.0, 0.0, 0.0);

        // Refine disk intersection point
        p_hit.z = self.height;

        // Compute error bounds for disk intersection
        let p_error = vector3(0.0, 0.0, 0.0);

        // Initialize SurfaceInteraction from parametric information
        let si = surface_interaction(
            p_hit,
            p_error,
            point2(u, v),
            -ray.d,
            dpdu,
            dpdv,
            dndu,
            dndv,
            ray.time,
            Some(Arc::new(self.clone())),
        );

        // Create hit.
        let isect = self.data.object_to_world.transform_surface_interaction(&si);
        Some(intersection(t_shape_hit, isect))
    }

    /// Returns `true` if a ray-shape intersection succeeds; otherwise `false`.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests (not supported).
    fn intersect_p(&self, r: &Ray, _test_alpha_texture: bool) -> bool {
        // Transform ray to object space
        //
        // We could just use transform_ray() but there is minor adjustment in
        // it that adjusts t_max which is not in transform_ray_with_error().
        let (ray, _o_err, _d_err) = self
            .data
            .world_to_object
            .clone()
            .unwrap()
            .transform_ray_with_error(r);

        // Compute plane intersection for disk.

        // Reject disk intersections for rays parallel to the disk's plane
        if ray.d.z == 0.0 {
            return false;
        }
        let t_shape_hit = (self.height - ray.o.z) / ray.d.z;
        if t_shape_hit <= 0.0 || t_shape_hit >= ray.t_max {
            return false;
        }

        // See if hit point is inside disk radii and phimax
        let p_hit = ray.at(t_shape_hit);
        let dist2 = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
        if dist2 > self.radius * self.radius || dist2 < self.inner_radius * self.inner_radius {
            return false;
        }

        // Test disk phi value against phimax
        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += TWO_PI;
        }
        if phi > self.phi_max {
            return false;
        }

        true
    }

    /// Returns the surface area of the shape in object space.
    fn area(&self) -> Float {
        self.phi_max * 0.5 * (self.radius * self.radius - self.inner_radius * self.inner_radius)
    }
}
