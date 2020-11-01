//! Hyperboloids

#![allow(dead_code)]
use super::{
    clamp, max, min, ArcTransform, Bounds3f, Dot, EFloat, Float, Intersection, Normal3, Point2f,
    Point3, Point3f, Quadratic, Ray, Shape, ShapeData, SurfaceInteraction, Vector3, TWO_PI,
};
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
    pub data: ShapeData,

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
pub fn hyperboloid(
    object_to_world: ArcTransform,
    world_to_object: ArcTransform,
    reverse_orientation: bool,
    point1: Point3f,
    point2: Point3f,
    phi_max: Float,
) -> Hyperboloid {
    let mut p1 = point1;
    let mut p2 = point2;

    let radius1 = (p1.x * p1.x + p1.y * p1.y).sqrt();
    let radius2 = (p2.x * p2.x + p2.y * p2.y).sqrt();

    let r_max = max(radius1, radius2);
    let z_min = min(p1.z, p2.z);
    let z_max = max(p1.z, p2.z);

    // Compute implicit function coefficients for hyperboloid
    if p2.z == 0.0 {
        swap(&mut p1, &mut p2);
    }

    let mut pp = p1;
    let mut xy1: Float;
    let mut xy2: Float;
    let mut ah: Float;
    let mut ch: Float;
    loop {
        pp += 2.0 * (p2 - p1);
        xy1 = pp.x * pp.x + pp.y * pp.y;
        xy2 = p2.x * p2.x + p2.y * p2.y;
        ah = (1.0 / xy1 - (pp.z * pp.z) / (xy1 * p2.z * p2.z))
            / (1.0 - (xy2 * pp.z * pp.z) / (xy1 * p2.z * p2.z));
        ch = (ah * xy2 - 1.0) / (p2.z * p2.z);

        if ah.is_finite() || ah != Float::NAN {
            break;
        }
    }

    Hyperboloid {
        p1,
        p2,
        z_min,
        z_max,
        r_max,
        ah,
        ch,
        phi_max: clamp(phi_max, 0.0, 360.0).to_radians(),
        data: ShapeData::new(
            object_to_world.clone(),
            Some(world_to_object.clone()),
            reverse_orientation,
        ),
    }
}

impl Shape for Hyperboloid {
    /// Returns the underlying shape data.
    fn get_data(&self) -> ShapeData {
        self.data.clone()
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
    fn intersect<'a>(&self, r: &Ray, _test_alpha_texture: bool) -> Option<Intersection<'a>> {
        // Transform ray to object space
        let (ray, o_err, d_err) = self
            .data
            .world_to_object
            .clone()
            .unwrap()
            .transform_ray_with_error(r);

        // Compute quadratic hyperboloid coefficients

        // Initialize EFloat ray coordinate values
        let ox = EFloat::new(ray.o.x, o_err.x);
        let oy = EFloat::new(ray.o.y, o_err.y);
        let oz = EFloat::new(ray.o.z, o_err.z);

        let dx = EFloat::new(ray.d.x, d_err.x);
        let dy = EFloat::new(ray.d.y, d_err.y);
        let dz = EFloat::new(ray.d.z, d_err.z);

        let a = self.ah * dx * dx + self.ah * dy * dy - self.ch * dz * dz;
        let b = 2.0 * (self.ah * dx * ox + self.ah * dy * oy - self.ch * dz * oz);
        let c = self.ah * ox * ox + self.ah * oy * oy - self.ch * oz * oz - 1.0;

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
                };
            }

            // Compute hyperboloid inverse mapping
            let mut p_hit = ray.at(Float::from(t_shape_hit));

            let mut v = (p_hit.z - self.p1.z) / (self.p2.z - self.p1.z);
            let mut pr = (1.0 - v) * self.p1 + v * self.p2;

            let mut phi = (pr.x * p_hit.y - p_hit.x * pr.y).atan2(p_hit.x * pr.x + p_hit.y * pr.y);
            if phi < 0.0 {
                phi += TWO_PI;
            }

            // Test hyperboloid intersection against clipping parameters
            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                if t_shape_hit == t1 {
                    return None;
                }

                t_shape_hit = t1;

                if t1.upper_bound() > ray.t_max {
                    return None;
                }

                // Compute hyperboloid inverse mapping
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

            // Find parametric representation of hyperboloid hit
            let u = phi / self.phi_max;

            // Compute hyperboloid dpdu and dpdv
            let cos_phi = phi.cos();
            let sin_phi = phi.sin();
            let dpdu = Vector3::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, 0.0);
            let dpdv = Vector3::new(
                (self.p2.x - self.p1.x) * cos_phi - (self.p2.y - self.p1.y) * sin_phi,
                (self.p2.x - self.p1.x) * sin_phi + (self.p2.y - self.p1.y) * cos_phi,
                self.p2.z - self.p1.z,
            );

            // Compute hyperboloid dndu and dndv
            let d2p_duu = -self.phi_max * self.phi_max * Vector3::new(p_hit.x, p_hit.y, 0.0);
            let d2p_duv = self.phi_max * Vector3::new(-dpdv.y, dpdv.x, 0.0);
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

            // Compute error bounds for hyperboloid intersection
            // Compute error bounds for intersection computed with ray equation
            let px = ox + t_shape_hit * dx;
            let py = oy + t_shape_hit * dy;
            let pz = oz + t_shape_hit * dz;
            let p_error = Vector3::new(
                px.get_absolute_error(),
                py.get_absolute_error(),
                pz.get_absolute_error(),
            );

            // Initialize SurfaceInteraction from parametric information
            let si = SurfaceInteraction::new(
                p_hit,
                p_error,
                Point2f::new(u, v),
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

        // Compute quadratic hyperboloid coefficients

        // Initialize EFloat ray coordinate values
        let ox = EFloat::new(ray.o.x, o_err.x);
        let oy = EFloat::new(ray.o.y, o_err.y);
        let oz = EFloat::new(ray.o.z, o_err.z);

        let dx = EFloat::new(ray.d.x, d_err.x);
        let dy = EFloat::new(ray.d.y, d_err.y);
        let dz = EFloat::new(ray.d.z, d_err.z);

        let a = self.ah * dx * dx + self.ah * dy * dy - self.ch * dz * dz;
        let b = 2.0 * (self.ah * dx * ox + self.ah * dy * oy - self.ch * dz * oz);
        let c = self.ah * ox * ox + self.ah * oy * oy - self.ch * oz * oz - 1.0;

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
                };
            }

            // Compute hyperboloid inverse mapping
            let mut p_hit = ray.at(Float::from(t_shape_hit));

            let mut v = (p_hit.z - self.p1.z) / (self.p2.z - self.p1.z);
            let mut pr = (1.0 - v) * self.p1 + v * self.p2;

            let mut phi = (pr.x * p_hit.y - p_hit.x * pr.y).atan2(p_hit.x * pr.x + p_hit.y * pr.y);
            if phi < 0.0 {
                phi += TWO_PI;
            }

            // Test hyperboloid intersection against clipping parameters
            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                if t_shape_hit == t1 {
                    return false;
                }

                t_shape_hit = t1;

                if t1.upper_bound() > ray.t_max {
                    return false;
                }

                // Compute hyperboloid inverse mapping
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
        } else {
            return false;
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
}
