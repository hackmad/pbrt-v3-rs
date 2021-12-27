//! Spheres

#![allow(dead_code)]
use core::efloat::*;
use core::geometry::*;
use core::interaction::*;
use core::paramset::*;
use core::pbrt::*;
use core::sampling::*;
use std::sync::Arc;

/// A sphere at origin [0, 0, 0].
#[derive(Clone)]
pub struct Sphere {
    /// Common shape data.
    pub data: Arc<ShapeData>,

    /// Radius of sphere.
    pub radius: Float,

    /// Minimum z-value to truncate sphere.
    pub z_min: Float,

    /// Maximum z-value to truncate sphere.
    pub z_max: Float,

    /// Minimum spherical coordinate for θ.
    pub theta_min: Float,

    /// Maximum spherical coordinate for θ.
    pub theta_max: Float,

    /// Maximum spherical coordinate for Φ.
    pub phi_max: Float,
}

impl Sphere {
    /// Create a new sphere at origin [0, 0, 0].
    ///
    /// * `object_to_world`     - The object to world transfomation.
    /// * `world_to_object`     - The world to object transfomation.
    /// * `reverse_orientation` - Indicates whether their surface normal directions
    ///                           should be reversed from the default
    /// * `radius`              - Radius of sphere.
    /// * `z_min`               - Minimum z-value to truncate sphere.
    /// * `z_max`               - Maximum z-value to truncate sphere.
    /// * `phi_max`             - Maximum spherical coordinate for Φ.
    pub fn new(
        object_to_world: ArcTransform,
        world_to_object: ArcTransform,
        reverse_orientation: bool,
        radius: Float,
        z_min: Float,
        z_max: Float,
        phi_max: Float,
    ) -> Self {
        Self {
            radius,
            z_min: clamp(min(z_min, z_max), -radius, radius),
            z_max: clamp(max(z_min, z_max), -radius, radius),
            theta_min: acos(clamp(min(z_min, z_max) / radius, -1.0, 1.0)),
            theta_max: acos(clamp(max(z_min, z_max) / radius, -1.0, 1.0)),
            phi_max: clamp(phi_max, 0.0, 360.0).to_radians(),
            data: Arc::new(ShapeData::new(
                Arc::clone(&object_to_world),
                Some(Arc::clone(&world_to_object)),
                reverse_orientation,
            )),
        }
    }
}

impl Shape for Sphere {
    /// Returns the shape type. Usually these are behind ArcShape and harder to
    /// debug. So this will be helpful.
    fn get_type(&self) -> &'static str {
        "sphere"
    }

    /// Returns the underlying shape data.
    fn get_data(&self) -> Arc<ShapeData> {
        Arc::clone(&self.data)
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
    fn intersect<'primitive, 'arena>(
        &self,
        r: &Ray,
        _test_alpha_texture: bool,
    ) -> Option<Intersection<'primitive, 'arena>> {
        // Transform ray to object space
        let (ray, o_err, d_err) = self
            .data
            .world_to_object
            .as_ref()
            .map(|w2o| w2o.transform_ray_with_error(r))
            .unwrap();

        // Compute quadratic sphere coefficients.

        // Initialize EFloat ray coordinate values.
        let ox = EFloat::new(ray.o.x, o_err.x);
        let oy = EFloat::new(ray.o.y, o_err.y);
        let oz = EFloat::new(ray.o.z, o_err.z);

        let dx = EFloat::new(ray.d.x, d_err.x);
        let dy = EFloat::new(ray.d.y, d_err.y);
        let dz = EFloat::new(ray.d.z, d_err.z);

        let a = dx * dx + dy * dy + dz * dz;
        let b = 2.0 * (dx * ox + dy * oy + dz * oz);
        let c = ox * ox + oy * oy + oz * oz - EFloat::from(self.radius) * EFloat::from(self.radius);

        // Solve quadratic equation for t values.
        let (t0, t1) = Quadratic::solve_efloat(a, b, c)?;

        // Check quadric shape t0 and t1 for nearest intersection.
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

        // Compute sphere hit position and phi.
        let mut p_hit = ray.at(Float::from(t_shape_hit));

        // Refine sphere intersection point.
        p_hit *= self.radius / p_hit.distance(Point3::new(0.0, 0.0, 0.0));

        if p_hit.x == 0.0 && p_hit.y == 0.0 {
            p_hit.x = 1e-5 * self.radius;
        }

        let mut phi = atan2(p_hit.y, p_hit.x);
        if phi < 0.0 {
            phi += TWO_PI;
        }

        // Test sphere intersection against clipping parameters.
        if (self.z_min > -self.radius && p_hit.z < self.z_min)
            || (self.z_max < self.radius && p_hit.z > self.z_max)
            || phi > self.phi_max
        {
            if t_shape_hit == t1 {
                return None;
            }
            if t1.upper_bound() > ray.t_max {
                return None;
            }

            t_shape_hit = t1;

            // Compute sphere hit position and phi.
            p_hit = ray.at(Float::from(t_shape_hit));

            // Refine sphere intersection point.
            p_hit *= self.radius / p_hit.distance(Point3::new(0.0, 0.0, 0.0));
            if p_hit.x == 0.0 && p_hit.y == 0.0 {
                p_hit.x = 1e-5 * self.radius;
            }

            phi = atan2(p_hit.y, p_hit.x);
            if phi < 0.0 {
                phi += TWO_PI;
            }

            if (self.z_min > -self.radius && p_hit.z < self.z_min)
                || (self.z_max < self.radius && p_hit.z > self.z_max)
                || phi > self.phi_max
            {
                return None;
            }
        }

        // Find parametric representation of sphere hit.
        let u = phi / self.phi_max;
        let theta = acos(clamp(p_hit.z / self.radius, -1.0, 1.0));
        let v = (theta - self.theta_min) / (self.theta_max - self.theta_min);

        // Compute sphere dpdu and dpdv.
        let z_radius = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
        let inv_z_radius = 1.0 / z_radius;
        let cos_phi = p_hit.x * inv_z_radius;
        let sin_phi = p_hit.y * inv_z_radius;
        let dpdu = Vector3::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, 0.0);
        let dpdv = (self.theta_max - self.theta_min)
            * Vector3::new(
                p_hit.z * cos_phi,
                p_hit.z * sin_phi,
                -self.radius * theta.sin(),
            );

        // Compute sphere dndu and dndv
        let d2p_duu = -self.phi_max * self.phi_max * Vector3::new(p_hit.x, p_hit.y, 0.0);
        let d2p_duv = (self.theta_max - self.theta_min)
            * p_hit.z
            * self.phi_max
            * Vector3::new(-sin_phi, cos_phi, 0.0);
        let d2p_dvv = -(self.theta_max - self.theta_min)
            * (self.theta_max - self.theta_min)
            * Vector3::new(p_hit.x, p_hit.y, p_hit.z);

        // Compute normal.
        let n = dpdu.cross(&dpdv).normalize();

        // Compute coefficients for first fundamental form.
        let e1 = dpdu.dot(&dpdu);
        let f1 = dpdu.dot(&dpdv);
        let g1 = dpdv.dot(&dpdv);

        // Compute coefficients for second fundamental form.
        let e2 = n.dot(&d2p_duu);
        let f2 = n.dot(&d2p_duv);
        let g2 = n.dot(&d2p_dvv);

        // Compute dndu and dndv from fundamental form coefficients.
        let inv_egf_1 = 1.0 / (e1 * g1 - f1 * f1);
        let dndu = Normal3::from(
            (f2 * f1 - e2 * g1) * inv_egf_1 * dpdu + (e2 * f1 - f2 * e1) * inv_egf_1 * dpdv,
        );
        let dndv = Normal3::from(
            (g2 * f1 - f2 * g1) * inv_egf_1 * dpdu + (f2 * f1 - g2 * e1) * inv_egf_1 * dpdv,
        );

        // Compute error bounds for sphere intersection
        let p_error = gamma(5) * Vector3::from(p_hit).abs();

        // Initialize SurfaceInteraction from parametric information.
        let mut si = SurfaceInteraction::new(
            p_hit,
            p_error,
            Point2::new(u, v),
            -ray.d,
            dpdu,
            dpdv,
            dndu,
            dndv,
            ray.time,
            Arc::clone(&self.data),
            0,
        );

        // Create hit.
        self.data
            .object_to_world
            .transform_surface_interaction(&mut si);

        Some(Intersection::new(Float::from(t_shape_hit), si))
    }

    /// Returns `true` if a ray-shape intersection succeeds; otherwise `false`.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests (not supported).
    fn intersect_p(&self, r: &Ray, _test_alpha_texture: bool) -> bool {
        // Transform ray to object space.
        let (ray, o_err, d_err) = self
            .data
            .world_to_object
            .as_ref()
            .map(|w2o| w2o.transform_ray_with_error(r))
            .unwrap();

        // Compute quadratic sphere coefficients.

        // Initialize EFloat ray coordinate values.
        let ox = EFloat::new(ray.o.x, o_err.x);
        let oy = EFloat::new(ray.o.y, o_err.y);
        let oz = EFloat::new(ray.o.z, o_err.z);

        let dx = EFloat::new(ray.d.x, d_err.x);
        let dy = EFloat::new(ray.d.y, d_err.y);
        let dz = EFloat::new(ray.d.z, d_err.z);

        let a = dx * dx + dy * dy + dz * dz;
        let b = 2.0 * (dx * ox + dy * oy + dz * oz);
        let c = ox * ox + oy * oy + oz * oz - EFloat::from(self.radius) * EFloat::from(self.radius);

        // Solve quadratic equation for `t` values.
        let solution = Quadratic::solve_efloat(a, b, c);
        if solution.is_none() {
            return false;
        }
        let (t0, t1) = solution.unwrap();

        // Check quadric shape _t0_ and _t1_ for nearest intersection
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

        // Compute sphere hit position and phi.
        let mut p_hit = ray.at(Float::from(t_shape_hit));

        // Refine sphere intersection point.
        p_hit *= self.radius / p_hit.distance(Point3::new(0.0, 0.0, 0.0));

        if p_hit.x == 0.0 && p_hit.y == 0.0 {
            p_hit.x = 1e-5 * self.radius;
        }

        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += TWO_PI;
        }

        // Test sphere intersection against clipping parameters.
        if (self.z_min > -self.radius && p_hit.z < self.z_min)
            || (self.z_max < self.radius && p_hit.z > self.z_max)
            || phi > self.phi_max
        {
            if t_shape_hit == t1 {
                return false;
            }
            if t1.upper_bound() > ray.t_max {
                return false;
            }

            t_shape_hit = t1;

            // Compute sphere hit position and phi.
            p_hit = ray.at(Float::from(t_shape_hit));

            // Refine sphere intersection point.
            p_hit *= self.radius / p_hit.distance(Point3::new(0.0, 0.0, 0.0));
            if p_hit.x == 0.0 && p_hit.y == 0.0 {
                p_hit.x = 1e-5 * self.radius;
            }

            phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 {
                phi += TWO_PI;
            }

            if (self.z_min > -self.radius && p_hit.z < self.z_min)
                || (self.z_max < self.radius && p_hit.z > self.z_max)
                || phi > self.phi_max
            {
                return false;
            }
        }

        true
    }

    /// Returns the surface area of the shape in object space.
    fn area(&self) -> Float {
        self.phi_max * self.radius * (self.z_max - self.z_min)
    }

    /// Sample a point on the surface and return the PDF with respect to area on
    /// the surface.
    ///
    /// NOTE: The returned `Hit` value will have `wo` = Vector3f::default().
    ///
    /// * `u` - Sample value to use.
    fn sample_area(&self, u: &Point2f) -> (Hit, Float) {
        let mut p_obj = Point3f::default() + self.radius * uniform_sample_sphere(u);

        let mut n = self
            .data
            .object_to_world
            .transform_normal(&Normal3f::new(p_obj.x, p_obj.y, p_obj.z))
            .normalize();
        if self.data.reverse_orientation {
            n *= -1.0;
        }

        // Reproject `p_obj` to sphere surface and compute `p_obj_error`.
        p_obj *= self.radius / p_obj.distance(Point3f::default());
        let p_obj_error = gamma(5) * Vector3f::from(p_obj).abs();
        let (p, p_error) = self
            .data
            .object_to_world
            .transform_point_with_abs_error(&p_obj, &p_obj_error);
        let it = Hit::new(p, 0.0, p_error, Vector3f::default(), n, None);
        let pdf = 1.0 / self.area();
        (it, pdf)
    }
}

impl From<(&ParamSet, ArcTransform, ArcTransform, bool)> for Sphere {
    /// Create a `Sphere` from given parameter set, object to world transform,
    /// world to object transform and whether or not surface normal orientation
    /// is reversed.
    ///
    /// * `p` - A tuple containing the parameter set, object to world transform,
    ///         world to object transform and whether or not surface normal
    ///         orientation is reversed.
    fn from(p: (&ParamSet, ArcTransform, ArcTransform, bool)) -> Self {
        let (params, o2w, w2o, reverse_orientation) = p;

        let radius = params.find_one_float("radius", 1.0);
        let z_min = params.find_one_float("zmin", -radius);
        let z_max = params.find_one_float("zmax", radius);
        let phi_max = params.find_one_float("phimax", 360.0);

        Self::new(
            Arc::clone(&o2w),
            Arc::clone(&w2o),
            reverse_orientation,
            radius,
            z_min,
            z_max,
            phi_max,
        )
    }
}
