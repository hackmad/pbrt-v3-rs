//! Hyperboloids

#![allow(dead_code)]
use core::efloat::*;
use core::geometry::*;
use core::interaction::*;
use core::paramset::*;
use core::pbrt::*;
use std::mem::swap;
use std::sync::Arc;

macro_rules! sqr {
    ($a: expr) => {
        $a * $a
    };
}

macro_rules! quad {
    ($a: expr) => {
        sqr!($a) * sqr!($a)
    };
}

/// A hyperboloid centered on the z-axis defined as a revolution of a line segment.
#[derive(Clone)]
pub struct Hyperboloid {
    /// Common shape data.
    pub data: Arc<ShapeData>,

    /// First point defining line segment to revolve.
    pub p1: Point3f,

    /// Second point defining line segment to revolve.
    pub p2: Point3f,

    /// Minimum z-value to truncate hyperboloid.
    pub z_min: Float,

    /// Maximum z-value to truncate hyperboloid.
    pub z_max: Float,

    /// Maximum angle Φ to truncate hyperboloid.
    pub phi_max: Float,

    /// Maximum of the 2 radii around z_min and z_max.
    pub r_max: Float,

    /// Implicit function coefficient.
    pub ah: Float,

    /// Implicit function coefficient.
    pub ch: Float,
}

impl Hyperboloid {
    /// Create a new hyperboloid centered on the z-axis defined as a revolution of
    /// a line segment.
    ///
    /// * `object_to_world`     - The object to world transfomation.
    /// * `world_to_object`     - The world to object transfomation.
    /// * `reverse_orientation` - Indicates whether their surface normal directions
    ///                           should be reversed from the default
    /// * `point1`              - First point defining line segment to revolve.
    /// * `point2`              - Second point defining line segment to revolve.
    /// * `phi_max`             - Maximum spherical coordinate for Φ.
    pub fn new(
        object_to_world: ArcTransform,
        world_to_object: ArcTransform,
        reverse_orientation: bool,
        point1: Point3f,
        point2: Point3f,
        phi_max: Float,
    ) -> Self {
        let mut p1 = point1;
        let mut p2 = point2;

        let radius1 = (p1.x * p1.x + p1.y * p1.y).sqrt();
        let radius2 = (p2.x * p2.x + p2.y * p2.y).sqrt();

        let r_max = max(radius1, radius2);
        let z_min = min(p1.z, p2.z);
        let z_max = max(p1.z, p2.z);

        // Compute implicit function coefficients for hyperboloid.
        if p2.z == 0.0 {
            swap(&mut p1, &mut p2);
        }

        let mut pp = p1;
        let mut ah: Float;
        let mut ch: Float;
        let mut count = 0;
        loop {
            pp += 2.0 * (p2 - p1);
            let xy1 = pp.x * pp.x + pp.y * pp.y;
            let xy2 = p2.x * p2.x + p2.y * p2.y;
            ah = (1.0 / xy1 - (pp.z * pp.z) / (xy1 * p2.z * p2.z))
                / (1.0 - (xy2 * pp.z * pp.z) / (xy1 * p2.z * p2.z));
            ch = (ah * xy2 - 1.0) / (p2.z * p2.z);

            if ah.is_finite() && !ah.is_nan() {
                break;
            }

            // Give up after 100,000 interations for fear of getting stuck in
            // infinte loop :(
            count += 1;
            if count > 100000 {
                panic!("Hyperboloid parameters possibly causing infinte loop.");
            }
        }

        Self {
            p1,
            p2,
            z_min,
            z_max,
            r_max,
            ah,
            ch,
            phi_max: clamp(phi_max, 0.0, 360.0).to_radians(),
            data: Arc::new(ShapeData::new(
                Arc::clone(&object_to_world),
                Some(Arc::clone(&world_to_object)),
                reverse_orientation,
            )),
        }
    }
}

impl Shape for Hyperboloid {
    /// Returns the shape type. Usually these are behind ArcShape and harder to
    /// debug. So this will be helpful.
    fn get_type(&self) -> &'static str {
        "hyperboloid"
    }

    /// Returns the underlying shape data.
    fn get_data(&self) -> Arc<ShapeData> {
        Arc::clone(&self.data)
    }

    /// Returns a bounding box in the shapes object space.
    fn object_bound(&self) -> Bounds3f {
        Bounds3f::new(
            Point3::new(-self.r_max, -self.r_max, self.z_min),
            Point3::new(self.r_max, self.r_max, self.z_max),
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
    ) -> Option<Intersection<'scene>> {
        // Transform ray to object space.
        let (ray, o_err, d_err) = self
            .data
            .world_to_object
            .as_ref()
            .map(|w2o| w2o.transform_ray_with_error(r))
            .unwrap();

        // Compute quadratic hyperboloid coefficients.

        // Initialize EFloat ray coordinate values.
        let ox = EFloat::new(ray.o.x, o_err.x);
        let oy = EFloat::new(ray.o.y, o_err.y);
        let oz = EFloat::new(ray.o.z, o_err.z);

        let dx = EFloat::new(ray.d.x, d_err.x);
        let dy = EFloat::new(ray.d.y, d_err.y);
        let dz = EFloat::new(ray.d.z, d_err.z);

        let a = self.ah * dx * dx + self.ah * dy * dy - self.ch * dz * dz;
        let b = 2.0 * (self.ah * dx * ox + self.ah * dy * oy - self.ch * dz * oz);
        let c = self.ah * ox * ox + self.ah * oy * oy - self.ch * oz * oz - 1.0;

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

        // Compute hyperboloid inverse mapping.
        let mut p_hit = ray.at(Float::from(t_shape_hit));

        let mut v = (p_hit.z - self.p1.z) / (self.p2.z - self.p1.z);
        let mut pr = (1.0 - v) * self.p1 + v * self.p2;

        let mut phi = (pr.x * p_hit.y - p_hit.x * pr.y).atan2(p_hit.x * pr.x + p_hit.y * pr.y);
        if phi < 0.0 {
            phi += TWO_PI;
        }

        // Test hyperboloid intersection against clipping parameters.
        if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
            if t_shape_hit == t1 {
                return None;
            }

            t_shape_hit = t1;

            if t1.upper_bound() > ray.t_max {
                return None;
            }

            // Compute hyperboloid inverse mapping.
            p_hit = ray.at(Float::from(t_shape_hit));

            v = (p_hit.z - self.p1.z) / (self.p2.z - self.p1.z);
            pr = (1.0 - v) * self.p1 + v * self.p2;

            phi = (pr.x * p_hit.y - p_hit.x * pr.y).atan2(p_hit.x * pr.x + p_hit.y * pr.y);
            if phi < 0.0 {
                phi += TWO_PI;
            }

            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                return None;
            }
        }

        // Find parametric representation of hyperboloid hit.
        let u = phi / self.phi_max;

        // Compute hyperboloid dpdu and dpdv.
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();
        let dpdu = Vector3::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, 0.0);
        let dpdv = Vector3::new(
            (self.p2.x - self.p1.x) * cos_phi - (self.p2.y - self.p1.y) * sin_phi,
            (self.p2.x - self.p1.x) * sin_phi + (self.p2.y - self.p1.y) * cos_phi,
            self.p2.z - self.p1.z,
        );

        // Compute hyperboloid dndu and dndv.
        let d2p_duu = -self.phi_max * self.phi_max * Vector3::new(p_hit.x, p_hit.y, 0.0);
        let d2p_duv = self.phi_max * Vector3::new(-dpdv.y, dpdv.x, 0.0);
        let d2p_dvv = Vector3::new(0.0, 0.0, 0.0);

        // Compute normal
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

        // Compute error bounds for hyperboloid intersection.

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
            Point2f::new(u, v),
            -ray.d,
            dpdu,
            dpdv,
            dndu,
            dndv,
            ray.time,
            Arc::clone(&self.data),
            0,
        );
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

        // Compute quadratic hyperboloid coefficients.

        // Initialize EFloat ray coordinate values.
        let ox = EFloat::new(ray.o.x, o_err.x);
        let oy = EFloat::new(ray.o.y, o_err.y);
        let oz = EFloat::new(ray.o.z, o_err.z);

        let dx = EFloat::new(ray.d.x, d_err.x);
        let dy = EFloat::new(ray.d.y, d_err.y);
        let dz = EFloat::new(ray.d.z, d_err.z);

        let a = self.ah * dx * dx + self.ah * dy * dy - self.ch * dz * dz;
        let b = 2.0 * (self.ah * dx * ox + self.ah * dy * oy - self.ch * dz * oz);
        let c = self.ah * ox * ox + self.ah * oy * oy - self.ch * oz * oz - 1.0;

        // Solve quadratic equation for t values.
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

        // Compute hyperboloid inverse mapping.
        let mut p_hit = ray.at(Float::from(t_shape_hit));

        let mut v = (p_hit.z - self.p1.z) / (self.p2.z - self.p1.z);
        let mut pr = (1.0 - v) * self.p1 + v * self.p2;

        let mut phi = (pr.x * p_hit.y - p_hit.x * pr.y).atan2(p_hit.x * pr.x + p_hit.y * pr.y);
        if phi < 0.0 {
            phi += TWO_PI;
        }

        // Test hyperboloid intersection against clipping parameters.
        if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
            if t_shape_hit == t1 {
                return false;
            }

            t_shape_hit = t1;

            if t1.upper_bound() > ray.t_max {
                return false;
            }

            // Compute hyperboloid inverse mapping.
            p_hit = ray.at(Float::from(t_shape_hit));

            v = (p_hit.z - self.p1.z) / (self.p2.z - self.p1.z);
            pr = (1.0 - v) * self.p1 + v * self.p2;

            phi = (pr.x * p_hit.y - p_hit.x * pr.y).atan2(p_hit.x * pr.x + p_hit.y * pr.y);
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
        self.phi_max / 6.0
            * (2.0 * quad!(self.p1.x) - 2.0 * self.p1.x * self.p1.x * self.p1.x * self.p2.x
                + 2.0 * quad!(self.p2.x)
                + 2.0
                    * (self.p1.y * self.p1.y + self.p1.y * self.p2.y + self.p2.y * self.p2.y)
                    * (sqr!(self.p1.y - self.p2.y) + sqr!(self.p1.z - self.p2.z))
                + self.p2.x
                    * self.p2.x
                    * (5.0 * self.p1.y * self.p1.y + 2.0 * self.p1.y * self.p2.y
                        - 4.0 * self.p2.y * self.p2.y
                        + 2.0 * sqr!(self.p1.z - self.p2.z))
                + self.p1.x
                    * self.p1.x
                    * (-4.0 * self.p1.y * self.p1.y
                        + 2.0 * self.p1.y * self.p2.y
                        + 5.0 * self.p2.y * self.p2.y
                        + 2.0 * sqr!(self.p1.z - self.p2.z))
                - 2.0
                    * self.p1.x
                    * self.p2.x
                    * (self.p2.x * self.p2.x - self.p1.y * self.p1.y + 5.0 * self.p1.y * self.p2.y
                        - self.p2.y * self.p2.y
                        - self.p1.z * self.p1.z
                        + 2.0 * self.p1.z * self.p2.z
                        - self.p2.z * self.p2.z))
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

impl From<(&ParamSet, ArcTransform, ArcTransform, bool)> for Hyperboloid {
    /// Create a `Hyperboloid` from given parameter set, object to world transform,
    /// world to object transform and whether or not surface normal orientation
    /// is reversed.
    ///
    /// * `p` - A tuple containing the parameter set, object to world transform,
    ///         world to object transform and whether or not surface normal
    ///         orientation is reversed.
    fn from(p: (&ParamSet, ArcTransform, ArcTransform, bool)) -> Self {
        let (params, o2w, w2o, reverse_orientation) = p;

        // NOTE: These are the original defaults in PBRT v2/v3 that cause
        // the constructor to get stuck in an infinite loop.
        //let p1 = params.find_one_point3f("p1", Point3f::new(0.0, 0.0, 0.0));
        //let p2 = params.find_one_point3f("p2", Point3f::new(1.0, 1.0, 1.0));

        let p1 = params.find_one_point3f("p1", Point3f::new(3.0, 1.2, 1.0));
        let p2 = params.find_one_point3f("p2", Point3f::new(1.0, -0.5, -1.0));

        let phi_max = params.find_one_float("phimax", 360.0);

        Self::new(
            Arc::clone(&o2w),
            Arc::clone(&w2o),
            reverse_orientation,
            p1,
            p2,
            phi_max,
        )
    }
}
