//! Paraboloids

#![allow(dead_code)]
use core::efloat::*;
use core::geometry::*;
use core::interaction::*;
use core::paramset::*;
use core::pbrt::*;
use std::sync::Arc;

/// A paraboloid centered on the z-axis.
#[derive(Clone)]
pub struct Paraboloid {
    /// Common shape data.
    pub data: Arc<ShapeData>,

    /// Radius of paraboloid.
    pub radius: Float,

    /// Minimum z-value to truncate paraboloid.
    pub z_min: Float,

    /// Maximum z-value to truncate paraboloid.
    pub z_max: Float,

    /// Maximum angle Φ to truncate paraboloid.
    pub phi_max: Float,
}

impl Paraboloid {
    /// Create a new paraboloid centered on the z-axis and base at origin [0, 0, 0].
    ///
    /// * `object_to_world`     - The object to world transfomation.
    /// * `world_to_object`     - The world to object transfomation.
    /// * `reverse_orientation` - Indicates whether their surface normal directions
    ///                           should be reversed from the default
    /// * `radius`              - Radius of paraboloid.
    /// * `z_min`               - Minimum z-value to truncate paraboloid.
    /// * `z_max`               - Maximum z-value to truncate paraboloid.
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
            z_min: min(z_min, z_max),
            z_max: max(z_min, z_max),
            phi_max: clamp(phi_max, 0.0, 360.0).to_radians(),
            data: Arc::new(ShapeData::new(
                Arc::clone(&object_to_world),
                Some(Arc::clone(&world_to_object)),
                reverse_orientation,
            )),
        }
    }
}

impl Shape for Paraboloid {
    /// Returns the shape type. Usually these are behind ArcShape and harder to
    /// debug. So this will be helpful.
    fn get_type(&self) -> &'static str {
        "paraboloid"
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
    fn intersect<'scene, 'arena>(
        &self,
        r: &Ray,
        _test_alpha_texture: bool,
    ) -> Option<Intersection<'scene, 'arena>> {
        // Transform ray to object space
        let (ray, o_err, d_err) = self
            .data
            .world_to_object
            .as_ref()
            .map(|w2o| w2o.transform_ray_with_error(r))
            .unwrap();

        // Compute quadratic paraboloid coefficients.

        // Initialize EFloat ray coordinate values.
        let ox = EFloat::new(ray.o.x, o_err.x);
        let oy = EFloat::new(ray.o.y, o_err.y);
        let oz = EFloat::new(ray.o.z, o_err.z);

        let dx = EFloat::new(ray.d.x, d_err.x);
        let dy = EFloat::new(ray.d.y, d_err.y);
        let dz = EFloat::new(ray.d.z, d_err.z);

        let k = EFloat::from(self.z_max) / (EFloat::from(self.radius) * EFloat::from(self.radius));

        let a = k * (dx * dx + dy * dy);
        let b = 2.0 * k * (dx * ox + dy * oy) - dz;
        let c = k * (ox * ox + oy * oy) - oz;

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
            };
        }

        // Compute paraboloid inverse mapping.
        let mut p_hit = ray.at(Float::from(t_shape_hit));

        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += TWO_PI;
        }

        // Test paraboloid intersection against clipping parameters.
        if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
            if t_shape_hit == t1 {
                return None;
            }

            t_shape_hit = t1;

            if t1.upper_bound() > ray.t_max {
                return None;
            }

            // Compute paraboloid inverse mapping.
            p_hit = ray.at(Float::from(t_shape_hit));

            phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 {
                phi += TWO_PI;
            }

            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                return None;
            }
        }

        // Find parametric representation of paraboloid hit.
        let u = phi / self.phi_max;
        let v = (p_hit.z - self.z_min) / (self.z_max - self.z_min);

        // Compute paraboloid dpdu and dpdv.
        let dpdu = Vector3::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, 0.0);
        let dpdv = (self.z_max - self.z_min)
            * Vector3::new(p_hit.x / (2.0 * p_hit.z), p_hit.y / (2.0 * p_hit.z), 1.0);

        // Compute paraboloid dndu and dndv.
        let d2p_duu = -self.phi_max * self.phi_max * Vector3::new(p_hit.x, p_hit.y, 0.0);
        let d2p_duv = (self.z_max - self.z_min)
            * self.phi_max
            * Vector3::new(-p_hit.y / (2.0 * p_hit.z), p_hit.x / (2.0 * p_hit.z), 0.0);
        let d2p_dvv = -(self.z_max - self.z_min)
            * (self.z_max - self.z_min)
            * Vector3::new(
                p_hit.x / (4.0 * p_hit.z * p_hit.z),
                p_hit.y / (4.0 * p_hit.z * p_hit.z),
                0.0,
            );

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
        let inv_e1g1f1_2 = 1.0 / (e1 * g1 - f1 * f1);
        let dndu = Normal3::from(
            (f2 * f1 - e2 * g1) * inv_e1g1f1_2 * dpdu + (e2 * f1 - f2 * e1) * inv_e1g1f1_2 * dpdv,
        );
        let dndv = Normal3::from(
            (g2 * f1 - f2 * g1) * inv_e1g1f1_2 * dpdu + (f2 * f1 - g2 * e1) * inv_e1g1f1_2 * dpdv,
        );

        // Compute error bounds for paraboloid intersection.

        // Compute error bounds for intersection computed with ray equation.
        let px = ox + t_shape_hit * dx;
        let py = oy + t_shape_hit * dy;
        let pz = oz + t_shape_hit * dz;
        let p_error = Vector3::new(
            px.get_absolute_error(),
            py.get_absolute_error(),
            pz.get_absolute_error(),
        );

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
        // Transform ray to object space
        let (ray, o_err, d_err) = self
            .data
            .world_to_object
            .as_ref()
            .map(|w2o| w2o.transform_ray_with_error(r))
            .unwrap();

        // Compute quadratic paraboloid coefficients.

        // Initialize EFloat ray coordinate values.
        let ox = EFloat::new(ray.o.x, o_err.x);
        let oy = EFloat::new(ray.o.y, o_err.y);
        let oz = EFloat::new(ray.o.z, o_err.z);

        let dx = EFloat::new(ray.d.x, d_err.x);
        let dy = EFloat::new(ray.d.y, d_err.y);
        let dz = EFloat::new(ray.d.z, d_err.z);

        let k = EFloat::from(self.z_max) / (EFloat::from(self.radius) * EFloat::from(self.radius));

        let a = k * (dx * dx + dy * dy);
        let b = 2.0 * k * (dx * ox + dy * oy) - dz;
        let c = k * (ox * ox + oy * oy) - oz;

        // Solve quadratic equation for t values
        let solution = Quadratic::solve_efloat(a, b, c);
        if solution.is_none() {
            return false;
        }
        let (t0, t1) = solution.unwrap();

        // Check quadric shape t0 and t1 for nearest intersection.
        if t0.upper_bound() > ray.t_max || t1.lower_bound() <= 0.0 {
            return false;
        }

        let mut t_shape_hit = t0;
        if t_shape_hit.lower_bound() <= 0.0 {
            t_shape_hit = t1;
            if t_shape_hit.upper_bound() > ray.t_max {
                return false;
            };
        }

        // Compute paraboloid inverse mapping.
        let mut p_hit = ray.at(Float::from(t_shape_hit));

        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += TWO_PI;
        }

        // Test paraboloid intersection against clipping parameters.
        if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
            if t_shape_hit == t1 {
                return false;
            }

            t_shape_hit = t1;

            if t1.upper_bound() > ray.t_max {
                return false;
            }

            // Compute paraboloid inverse mapping.
            p_hit = ray.at(Float::from(t_shape_hit));

            phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 {
                phi += TWO_PI;
            }

            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                return false;
            }
        }

        true
    }

    /// Returns the surface area of the shape in object space.
    fn area(&self) -> Float {
        let radius2 = self.radius * self.radius;
        let k = 4.0 * self.z_max / radius2;
        (radius2 * radius2 * self.phi_max / (12.0 * self.z_max * self.z_max))
            * ((k * self.z_max + 1.0).powf(1.5) - (k * self.z_min + 1.0).powf(1.5))
    }

    /// Sample a point on the surface and return the PDF with respect to area on
    /// the surface.
    ///
    /// NOTE: The returned `Hit` value will have `wo` = Vector3f::ZERO.
    ///
    /// * `u` - Sample value to use.
    fn sample_area(&self, _u: &Point2f) -> (Hit, Float) {
        todo!()
    }
}

impl From<(&ParamSet, ArcTransform, ArcTransform, bool)> for Paraboloid {
    /// Create a `Paraboloid` from given parameter set, object to world transform,
    /// world to object transform and whether or not surface normal orientation
    /// is reversed.
    ///
    /// * `p` - A tuple containing the parameter set, object to world transform,
    ///         world to object transform and whether or not surface normal
    ///         orientation is reversed.
    fn from(p: (&ParamSet, ArcTransform, ArcTransform, bool)) -> Self {
        let (params, o2w, w2o, reverse_orientation) = p;

        let radius = params.find_one_float("radius", 1.0);
        let z_min = params.find_one_float("zmin", 0.0);
        let z_max = params.find_one_float("zmax", 1.0);
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
