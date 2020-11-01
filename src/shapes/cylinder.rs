//! Cylinders

#![allow(dead_code)]
use super::{
    clamp, gamma, max, min, ArcTransform, Bounds3f, Dot, EFloat, Float, Intersection, Normal3,
    Point2, Point3, Quadratic, Ray, Shape, ShapeData, SurfaceInteraction, Vector3, TWO_PI,
};
use std::sync::Arc;

/// A cylinder.
#[derive(Clone)]
pub struct Cylinder {
    /// Common shape data.
    pub data: ShapeData,

    /// Radius of cylinder.
    pub radius: Float,

    /// Minimum z-value to truncate cylinder.
    pub z_min: Float,

    /// Maximum z-value to truncate cylinder.
    pub z_max: Float,

    /// Maximum angle Φ to truncate cylinder.
    pub phi_max: Float,
}

/// Create a new cylinder.
///
/// * `object_to_world`     - The object to world transfomation.
/// * `world_to_object`     - The world to object transfomation.
/// * `reverse_orientation` - Indicates whether their surface normal directions
///                           should be reversed from the default
/// * `radius`              - Radius of cylinder.
/// * `z_min`               - Minimum z-value to truncate cylinder.
/// * `z_max`               - Maximum z-value to truncate cylinder.
/// * `phi_max`             - Maximum angle Φ to truncate cylinder.
pub fn cylinder(
    object_to_world: ArcTransform,
    world_to_object: ArcTransform,
    reverse_orientation: bool,
    radius: Float,
    z_min: Float,
    z_max: Float,
    phi_max: Float,
) -> Cylinder {
    let zmin = clamp(min(z_min, z_max), -radius, radius);
    let zmax = clamp(max(z_min, z_max), -radius, radius);
    Cylinder {
        radius,
        z_min: zmin,
        z_max: zmax,
        phi_max: clamp(phi_max, 0.0, 360.0).to_radians(),
        data: ShapeData::new(
            object_to_world.clone(),
            Some(world_to_object.clone()),
            reverse_orientation,
        ),
    }
}

impl Shape for Cylinder {
    /// Returns the underlying shape data.
    fn get_data(&self) -> ShapeData {
        self.data.clone()
    }

    /// Returns a bounding box in the shapes object space.
    fn object_bound(&self) -> Bounds3f {
        Bounds3f::new(
            Point3::new(-self.radius, -self.radius, self.z_min),
            Point3::new(self.radius, self.radius, self.z_max),
        )
    }

    /// Returns geometric details if a ray intersects the shape intersection.
    /// If there is no intersection, `None` is returned.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests (not supported).
    fn intersect<'a>(&self, r: &Ray, _test_alpha_texture: bool) -> Option<Intersection<'a>> {
        // Transform ray to object space
        let (ray, o_err, d_err) = self
            .data
            .world_to_object
            .clone()
            .unwrap()
            .transform_ray_with_error(r);

        // Compute quadratic cylinder coefficients

        // Initialize EFloat ray coordinate values
        let ox = EFloat::new(ray.o.x, o_err.x);
        let oy = EFloat::new(ray.o.y, o_err.y);

        let dx = EFloat::new(ray.d.x, d_err.x);
        let dy = EFloat::new(ray.d.y, d_err.y);

        let a = dx * dx + dy * dy;
        let b = 2.0 * (dx * ox + dy * oy);
        let c = ox * ox + oy * oy - EFloat::from(self.radius) * EFloat::from(self.radius);

        // Solve quadratic equation for t values
        if let Some((t0, t1)) = Quadratic::solve(a, b, c) {
            // Check quadric shape t0 and t1 for nearest intersection
            if t0.upper_bound() > ray.t_max || t1.lower_bound() <= 0.0 {
                return None;
            }

            let mut t_shape_hit = t0;
            if t_shape_hit.lower_bound() <= 0.0 {
                t_shape_hit = t1;
                if t_shape_hit.upper_bound() > ray.t_max {
                    return None;
                }
            }

            // Compute cylinder hit position and phi
            let mut p_hit = ray.at(Float::from(t_shape_hit));

            // Refine cylinder intersection point
            let hit_rad = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
            p_hit.x *= self.radius / hit_rad;
            p_hit.y *= self.radius / hit_rad;

            let mut phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 {
                phi += TWO_PI;
            }

            // Test cylinder intersection against clipping parameters
            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                if t_shape_hit == t1 {
                    return None;
                }
                if t1.upper_bound() > ray.t_max {
                    return None;
                }

                t_shape_hit = t1;

                // Compute cylinder hit position and phi
                p_hit = ray.at(Float::from(t_shape_hit));

                // Refine cylinder intersection point
                let hit_rad = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
                p_hit.x *= self.radius / hit_rad;
                p_hit.y *= self.radius / hit_rad;

                phi = p_hit.y.atan2(p_hit.x);
                if phi < 0.0 {
                    phi += TWO_PI;
                }

                if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                    return None;
                }
            }

            // Find parametric representation of cylinder hit
            let u = phi / self.phi_max;
            let v = (p_hit.z - self.z_min) / (self.z_max - self.z_min);

            // Compute cylinder dpdu and dpdv
            let dpdu = Vector3::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, 0.0);
            let dpdv = Vector3::new(0.0, 0.0, self.z_max - self.z_min);

            // Compute cylinder dndu and dndv
            let d2p_duu = -self.phi_max * self.phi_max * Vector3::new(p_hit.x, p_hit.y, 0.0);
            let d2p_duv = Vector3::new(0.0, 0.0, 0.0);
            let d2p_dvv = Vector3::new(0.0, 0.0, 0.0);

            // Compute normal
            let n = dpdu.cross(&dpdv).normalize();

            // Compute coefficients for first fundamental form
            let e1 = dpdu.dot(&dpdu);
            let f1 = dpdu.dot(&dpdv);
            let g1 = dpdv.dot(&dpdv);

            // Compute coefficients for second fundamental form.
            let e2 = n.dot(&d2p_duu);
            let f2 = n.dot(&d2p_duv);
            let g2 = n.dot(&d2p_dvv);

            // Compute dndu and dndv from fundamental form coefficients
            let inv_egf_1 = 1.0 / (e1 * g1 - f1 * f1);
            let dndu = Normal3::from(
                (f2 * f1 - e2 * g1) * inv_egf_1 * dpdu + (e2 * f1 - f2 * e1) * inv_egf_1 * dpdv,
            );
            let dndv = Normal3::from(
                (g2 * f1 - f2 * g1) * inv_egf_1 * dpdu + (f2 * f1 - g2 * e1) * inv_egf_1 * dpdv,
            );

            // Compute error bounds for cylinder intersection
            let p_error = gamma(3) * Vector3::new(p_hit.x, p_hit.y, 0.0).abs();

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
            let t_hit = Float::from(t_shape_hit);
            Some(Intersection::new(t_hit, isect))
        } else {
            None
        }
    }

    /// Returns `true` if a ray-shape intersection succeeds; otherwise `false`.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests (not supported).
    fn intersect_p(&self, r: &Ray, _test_alpha_texture: bool) -> bool {
        // Transform ray to object space
        let (ray, o_err, d_err) = self
            .data
            .world_to_object
            .clone()
            .unwrap()
            .transform_ray_with_error(r);

        // Compute quadratic cylinder coefficients

        // Initialize EFloat ray coordinate values
        let ox = EFloat::new(ray.o.x, o_err.x);
        let oy = EFloat::new(ray.o.y, o_err.y);

        let dx = EFloat::new(ray.d.x, d_err.x);
        let dy = EFloat::new(ray.d.y, d_err.y);

        let a = dx * dx + dy * dy;
        let b = 2.0 * (dx * ox + dy * oy);
        let c = ox * ox + oy * oy - EFloat::from(self.radius) * EFloat::from(self.radius);

        // Solve quadratic equation for t values
        if let Some((t0, t1)) = Quadratic::solve(a, b, c) {
            // Check quadric shape t0 and t1 for nearest intersection
            if t0.upper_bound() > ray.t_max || t1.lower_bound() <= 0.0 {
                return false;
            }

            let mut t_shape_hit = t0;
            if t_shape_hit.lower_bound() <= 0.0 {
                t_shape_hit = t1;
                if t_shape_hit.upper_bound() > ray.t_max {
                    return false;
                }
            }

            // Compute cylinder hit position and phi
            let mut p_hit = ray.at(Float::from(t_shape_hit));

            // Refine cylinder intersection point
            let hit_rad = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
            p_hit.x *= self.radius / hit_rad;
            p_hit.y *= self.radius / hit_rad;

            let mut phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 {
                phi += TWO_PI;
            }

            // Test cylinder intersection against clipping parameters
            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                if t_shape_hit == t1 {
                    return false;
                }
                if t1.upper_bound() > ray.t_max {
                    return false;
                }

                t_shape_hit = t1;

                // Compute cylinder hit position and phi
                p_hit = ray.at(Float::from(t_shape_hit));

                // Refine cylinder intersection point
                let hit_rad = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
                p_hit.x *= self.radius / hit_rad;
                p_hit.y *= self.radius / hit_rad;

                phi = p_hit.y.atan2(p_hit.x);
                if phi < 0.0 {
                    phi += TWO_PI;
                }

                if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                    return false;
                }
            }
        } else {
            return false;
        }

        true
    }

    /// Returns the surface area of the shape in object space.
    fn area(&self) -> Float {
        (self.z_max - self.z_min) * self.radius * self.phi_max
    }
}
