//! Disks

#![allow(dead_code)]
use crate::core::paramset::*;
use crate::core::pbrt::*;
use crate::core::{geometry::*, sampling::concentric_sample_disk};
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

impl Disk {
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
    pub fn new(
        object_to_world: ArcTransform,
        world_to_object: ArcTransform,
        reverse_orientation: bool,
        height: Float,
        radius: Float,
        inner_radius: Float,
        phi_max: Float,
    ) -> Self {
        Self {
            height,
            radius,
            inner_radius,
            phi_max: clamp(phi_max, 0.0, 360.0).to_radians(),
            data: ShapeData::new(
                object_to_world.clone(),
                Some(world_to_object.clone()),
                reverse_orientation,
            ),
        }
    }
}

impl Shape for Disk {
    /// Returns the underlying shape data.
    fn get_data(&self) -> ShapeData {
        self.data.clone()
    }

    /// Returns a bounding box in the shapes object space.
    fn object_bound(&self) -> Bounds3f {
        Bounds3f::new(
            Point3::new(-self.radius, -self.radius, self.height),
            Point3::new(self.radius, self.radius, self.height),
        )
    }

    /// Returns geometric details if a ray intersects the shape intersection.
    /// If there is no intersection, `None` is returned.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests (not supported).
    fn intersect<'a>(&self, r: &Ray, _test_alpha_texture: bool) -> Option<Intersection<'a>> {
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
        let dpdu = Vector3::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, 0.0);
        let dpdv = Vector3::new(p_hit.x, p_hit.y, 0.0) * (self.inner_radius - self.radius) / r_hit;
        let dndu = Normal3::new(0.0, 0.0, 0.0);
        let dndv = Normal3::new(0.0, 0.0, 0.0);

        // Refine disk intersection point
        p_hit.z = self.height;

        // Compute error bounds for disk intersection
        let p_error = Vector3::new(0.0, 0.0, 0.0);

        // Initialize SurfaceInteraction from parametric information
        let si = SurfaceInteraction::new(
            p_hit,
            p_error,
            Point2::new(u, v),
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
        Some(Intersection::new(t_shape_hit, isect))
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

    /// Sample a point on the surface and return the PDF with respect to area on
    /// the surface.
    ///
    /// NOTE: The returned `Hit` value will have `wo` = Vector3f::default().
    ///
    /// * `u` - Sample value to use.
    fn sample_area(&self, u: &Point2f) -> (Hit, Float) {
        let pd = concentric_sample_disk(u);
        let p_obj = Point3f::new(pd.x * self.radius, pd.y * self.radius, self.height);

        let mut n = self
            .data
            .object_to_world
            .transform_normal(&Normal3f::new(0.0, 0.0, 1.0))
            .normalize();
        if self.data.reverse_orientation {
            n *= -1.0;
        }

        let p = self.data.object_to_world.transform_point(&p_obj);
        let p_error = self
            .data
            .object_to_world
            .transform_point_abs_error(&p_obj, &Vector3f::default());
        let it = Hit::new(p, 0.0, p_error, Vector3f::default(), n, None);
        let pdf = 1.0 / self.area();
        (it, pdf)
    }

    /// Sample a point on the shape given a reference point and return the PDF
    /// with respect to the solid angle from ref.
    ///
    /// * `hit` - Reference point on shape.
    /// * `u`   - Sample value to use.
    fn sample_solid_angle(&self, hit: &Hit, u: &Point2f) -> (Hit, Float) {
        todo!()
    }

    /// Returns the PDF with respect to solid angle.
    ///
    /// * `hit` - The interaction hit point.
    /// * `wi`  - The incoming direction.
    fn pdf_solid_angle(&self, hit: &Hit, wi: &Vector3f) -> Float {
        todo!()
    }
}

impl From<(&ParamSet, ArcTransform, ArcTransform, bool)> for Disk {
    /// Create a `Disk` from given parameter set, object to world transform,
    /// world to object transform and whether or not surface normal orientation
    /// is reversed.
    ///
    /// * `p` - A tuple containing the parameter set, object to world transform,
    ///         world to object transform and whether or not surface normal
    ///         orientation is reversed.
    fn from(p: (&ParamSet, ArcTransform, ArcTransform, bool)) -> Self {
        let (params, o2w, w2o, reverse_orientation) = p;

        let height = params.find_one_float("height", 0.0);
        let radius = params.find_one_float("radius", 1.0);
        let inner_radius = params.find_one_float("innerradius", 0.0);
        let phi_max = params.find_one_float("phimax", 360.0);

        Self::new(
            o2w.clone(),
            w2o.clone(),
            reverse_orientation,
            height,
            radius,
            inner_radius,
            phi_max,
        )
    }
}
