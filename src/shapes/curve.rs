//! Curves

#![allow(dead_code)]
use super::{
    abs, bounds3, clamp, coordinate_system, intersection, lerp, log2, look_at, max, min, point2,
    rotate_axis, shape_data, surface_interaction, vector3, ArcShape, ArcTransform, Bounds3f, Dot,
    Float, Intersection, Normal3f, Point3f, Ray, Shape, ShapeData, Union, Vector3, Vector3f,
};
use std::sync::Arc;

/// Curve types.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum CurveType {
    /// Flat curve is always oriented perpendiculur to an approaching ray.
    Flat,

    /// Cylinder curve has normal shading and appears cylinderical.
    Cylinder,

    /// Ribbon curve has fixed orientation at start and end points;
    /// intermediate orientations are smoothly interpolated.
    Ribbon,
}

/// A curve modeled as a cubic Bézier spline given by the polynomial:
/// p(u) = (1 - u^3) * p0  + 3 * (1 -u^2) * u * p1 + 3 * (1 - u) * u^2 p2 + u^3 * p3
#[derive(Clone)]
pub struct Curve {
    /// Common shape data.
    pub data: ShapeData,

    /// Common curve parameters.
    pub common: CurveCommon,

    /// Minimum u-parameter for the curve.
    pub u_min: Float,

    /// Maximum u-parameter for the curve.
    pub u_max: Float,
}

/// Create a new curve.
///
/// * `object_to_world`     - The object to world transfomation.
/// * `world_to_object`     - The world to object transfomation.
/// * `reverse_orientation` - Indicates whether their surface normal directions
///                           should be reversed from the default
/// * `u_min`               - Minimum u-parameter for the curve.
/// * `u_max`               - Maximum u-parameter for the curve.
pub fn curve(
    object_to_world: ArcTransform,
    world_to_object: ArcTransform,
    reverse_orientation: bool,
    common: CurveCommon,
    u_min: Float,
    u_max: Float,
) -> Curve {
    Curve {
        data: shape_data(
            object_to_world.clone(),
            Some(world_to_object.clone()),
            reverse_orientation,
        ),
        common,
        u_min,
        u_max,
    }
}

/// Create curve segments.
///
/// * `object_to_world`     - The object to world transfomation.
/// * `world_to_object`     - The world to object transfomation.
/// * `reverse_orientation` - Indicates whether their surface normal directions
///                           should be reversed from the default
/// * `curve_type`          - Curve type.
/// * `c`                   - Object space control points.
/// * `width`               - The width of the curve at the start and end points.
/// * `norm`                - Surface normal at the start and end points.
/// * `split_depth`         - Split depth.
pub fn create_curve(
    object_to_world: ArcTransform,
    world_to_object: ArcTransform,
    reverse_orientation: bool,
    curve_type: CurveType,
    c: [Point3f; 4],
    width: [Float; 2],
    norm: [Normal3f; 2],
    split_depth: i32,
) -> Vec<ArcShape> {
    let common = curve_common(curve_type, c, width, norm);

    let num_segments = 1 << split_depth;
    let mut segments = Vec::<ArcShape>::with_capacity(num_segments);

    let f = 1.0 / num_segments as Float;
    for i in 0..num_segments {
        let u_min = i as Float * f;
        let u_max = (i + 1) as Float * f;
        let curve = curve(
            object_to_world.clone(),
            world_to_object.clone(),
            reverse_orientation,
            common,
            u_min,
            u_max,
        );
        segments.push(Arc::new(curve));
    }

    segments
}

impl Curve {
    /// Recursively split curve in 2 sections if there is an intersection to
    /// find the actual intersection.
    ///
    /// * `ray`           - The ray.
    /// * `cp`            - The control points.
    /// * `ray_to_object` - Transform to bring things out of the coordinate system
    ///                     centered at a ray's origin with the ray's direction
    ///                     as +z axis. Generally created with look_at transform.
    /// * `u0`            - The starting u-parameter.
    /// * `u1`            - The ending u-parameter.
    /// * `depth`         - The recursion depth.
    /// * `current_hit`   - The current hit found so far. Initialize to `None`
    ///                     for first call.
    fn recursive_intersect(
        &self,
        ray: &Ray,
        cp: &[Point3f],
        ray_to_object: ArcTransform,
        u0: Float,
        u1: Float,
        depth: u32,
        current_hit: Option<Intersection>,
    ) -> Option<Intersection> {
        let ray_length = ray.d.length();
        let z_max = ray_length * ray.t_max;

        if depth > 0 {
            // Split curve segment into sub-segments and test for intersection
            let cp_split = subdivide_bezier(cp);

            // For each of the two segments, see if the ray's bounding box
            // overlaps the segment before recursively checking for
            // intersection with it.
            let u = [u0, (u0 + u1) / 2.0, u1];
            for seg in 0..2 {
                // Splice containing the 4 control poitns for the current segment.
                let cps = &cp_split[seg * 3..seg * 3 + 4];

                let max_width = max(
                    lerp(u[seg], self.common.width[0], self.common.width[1]),
                    lerp(u[seg + 1], self.common.width[0], self.common.width[1]),
                );

                // As above, check y first, since it most commonly lets us exit
                // out early.
                if max(max(cps[0].y, cps[1].y), max(cps[2].y, cps[3].y)) + 0.5 * max_width < 0.0
                    || min(min(cps[0].y, cps[1].y), min(cps[2].y, cps[3].y)) - 0.5 * max_width > 0.0
                {
                    continue;
                }

                if max(max(cps[0].x, cps[1].x), max(cps[2].x, cps[3].x)) + 0.5 * max_width < 0.0
                    || min(min(cps[0].x, cps[1].x), min(cps[2].x, cps[3].x)) - 0.5 * max_width > 0.0
                {
                    continue;
                }

                if max(max(cps[0].z, cps[1].z), max(cps[2].z, cps[3].z)) + 0.5 * max_width < 0.0
                    || min(min(cps[0].z, cps[1].z), min(cps[2].z, cps[3].z)) - 0.5 * max_width
                        > z_max
                {
                    continue;
                }

                let hit = self.recursive_intersect(
                    ray,
                    cps,
                    ray_to_object.clone(),
                    u[seg],
                    u[seg + 1],
                    depth - 1,
                    current_hit.clone(),
                );

                // If we found an intersection we can exit out immediately.
                if !hit.is_none() {
                    return hit;
                }
            }

            return current_hit;
        }

        // Intersect ray with curve segment.

        // Test ray against segment endpoint boundaries.

        // Test sample point against tangent perpendicular at curve start.
        let mut edge = (cp[1].y - cp[0].y) * (-cp[0].y) + cp[0].x * (cp[0].x - cp[1].x);
        if edge < 0.0 {
            return None;
        }

        // Test sample point against tangent perpendicular at curve end.
        edge = (cp[2].y - cp[3].y) * (-cp[3].y) + cp[3].x * (cp[3].x - cp[2].x);
        if edge < 0.0 {
            return None;
        }

        // Compute line w that gives minimum distance to sample point.
        let segment_direction = cp[3] - cp[0];
        let denom = segment_direction.length_squared();
        if denom == 0.0 {
            return None;
        }
        let w = (-Vector3::from(cp[0])).dot(&segment_direction) / denom;

        // Compute u coordinate of curve intersection point and hit_width.
        let u = clamp(lerp(w, u0, u1), u0, u1);
        let mut hit_width = lerp(u, self.common.width[0], self.common.width[1]);
        let mut n_hit = Normal3f::default();
        if self.common.curve_type == CurveType::Ribbon {
            // Scale hit_width based on ribbon orientation.
            let sin0 =
                ((1.0 - u) * self.common.normal_angle).sin() * self.common.inv_sin_normal_angle;
            let sin1 = (u * self.common.normal_angle).sin() * self.common.inv_sin_normal_angle;
            n_hit = sin0 * self.common.n[0] + sin1 * self.common.n[1];
            hit_width *= n_hit.abs_dot(&ray.d) / ray_length;
        }

        // Test intersection point against curve width
        let (pc, dpcdw) = eval_bezier(cp, clamp(w, 0.0, 1.0));
        let pt_curve_dist2 = pc.x * pc.x + pc.y * pc.y;
        if pt_curve_dist2 > hit_width * hit_width * 0.25 {
            return None;
        }
        if pc.z < 0.0 || pc.z > z_max {
            return None;
        }

        // Compute v coordinate of curve intersection point
        let pt_curve_dist = pt_curve_dist2.sqrt();
        let edge_func = dpcdw.x * (-pc.y) + pc.x * dpcdw.y;
        let v = if edge_func > 0.0 {
            0.5 + pt_curve_dist / hit_width
        } else {
            0.5 - pt_curve_dist / hit_width
        };

        // Not sure if this check is even necessary.
        if current_hit.is_none() {
            // Compute hit `t` and partial derivatives for curve intersection.
            let t_hit = pc.z / ray_length;

            // Compute error bounds for curve intersection
            let p_error = vector3(2.0 * hit_width, 2.0 * hit_width, 2.0 * hit_width);

            // Compute dpdu and dpdv for curve intersection.
            let (_, dpdu) = eval_bezier(&self.common.cp_obj, u);
            assert!(
                dpdu != Vector3f::default(),
                "u={}, cp=[{:?}]",
                u,
                self.common.cp_obj
            );

            let dpdv = if self.common.curve_type == CurveType::Ribbon {
                Vector3::from(n_hit).cross(&dpdu).normalize() * hit_width
            } else {
                // Compute curve dpdv for flat and cylinder curves
                let dpdu_plane = ray_to_object.inverse().transform_vector(&dpdu);
                let mut dpdv_plane =
                    vector3(-dpdu_plane.y, dpdu_plane.x, 0.0).normalize() * hit_width;
                if self.common.curve_type == CurveType::Cylinder {
                    // Rotate dpdv_plane to give cylindrical appearance.
                    let theta = lerp(v, -90.0, 90.0);
                    let rot = rotate_axis(-theta, &dpdu_plane);
                    dpdv_plane = rot.transform_vector(&dpdv_plane);
                }
                ray_to_object.transform_vector(&dpdv_plane)
            };

            let si = surface_interaction(
                ray.at(t_hit),
                p_error,
                point2(u, v),
                -ray.d,
                dpdu,
                dpdv,
                Normal3f::default(),
                Normal3f::default(),
                ray.time,
                Some(Arc::new(self.clone())),
            );
            let isect = self
                .data
                .object_to_world
                .clone()
                .transform_surface_interaction(&si);

            Some(intersection(t_hit, isect))
        } else {
            current_hit
        }
    }

    /// Computes the blossom the spline.
    fn blossom_bezier(&self) -> [Point3f; 4] {
        [
            blossom_bezier(&self.common.cp_obj, self.u_min, self.u_min, self.u_min),
            blossom_bezier(&self.common.cp_obj, self.u_min, self.u_min, self.u_max),
            blossom_bezier(&self.common.cp_obj, self.u_min, self.u_max, self.u_max),
            blossom_bezier(&self.common.cp_obj, self.u_max, self.u_max, self.u_max),
        ]
    }
}

impl Shape for Curve {
    /// Returns the underlying shape data.
    fn get_data(&self) -> ShapeData {
        self.data.clone()
    }

    /// Returns a bounding box in the shapes object space.
    fn object_bound(&self) -> Bounds3f {
        // Compute object-space control points for curve segment, cp_obj
        let cp_obj = self.blossom_bezier();

        // Using the convex hull property; i.e. the curve lies within the convex
        // hull of its control points. Then expand the bounds by half the maximum
        // width over the entire parameteric extent of the curve.
        let width = [
            lerp(self.u_min, self.common.width[0], self.common.width[1]),
            lerp(self.u_max, self.common.width[0], self.common.width[1]),
        ];

        bounds3(cp_obj[0], cp_obj[1])
            .union(&bounds3(cp_obj[2], cp_obj[3]))
            .expand(max(width[0], width[1]) * 0.5)
    }

    /// Returns geometric details if a ray intersects the shape intersection.
    /// If there is no intersection, `None` is returned.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests (not supported).
    fn intersect(&self, r: &Ray, _test_alpha_texture: bool) -> Option<Intersection> {
        // Transform ray to object space.
        //
        // We could just use transform_ray() but there is minor adjustment in
        // it that adjusts t_max which is not in transform_ray_with_error().
        let (ray, _o_err, _d_err) = self
            .data
            .world_to_object
            .clone()
            .unwrap()
            .transform_ray_with_error(r);

        // Compute object-space control points for curve segment, cp_obj.
        let cp_obj = self.blossom_bezier();

        // Project curve control points to plane perpendicular to ray.

        // Be careful to set the "up" direction passed to LookAt() to equal the
        // vector from the first to the last control points.  In turn, this
        // helps orient the curve to be roughly parallel to the x axis in the
        // ray coordinate system.
        //
        // In turn (especially for curves that are approaching stright lines),
        // we get curve bounds with minimal extent in y, which in turn lets us
        // early out more quickly in recursiveIntersect().
        let mut dx = ray.d.cross(&(cp_obj[3] - cp_obj[0]));
        if dx.length_squared() == 0.0 {
            // If the ray and the vector between the first and last control
            // points are parallel, dx will be zero.  Generate an arbitrary xy
            // orientation for the ray coordinate system so that intersection
            // tests can proceeed in this unusual case.
            let mut dy = Vector3f::default();
            coordinate_system(&ray.d, &mut dx, &mut dy);
        }

        let object_to_ray = look_at(&ray.o, &(ray.o + ray.d), &dx);
        let cp = [
            object_to_ray.transform_point(&cp_obj[0]),
            object_to_ray.transform_point(&cp_obj[1]),
            object_to_ray.transform_point(&cp_obj[2]),
            object_to_ray.transform_point(&cp_obj[3]),
        ];

        // Before going any further, see if the ray's bounding box intersects
        // the curve's bounding box. We start with the y dimension, since the y
        // extent is generally the smallest (and is often tiny) due to our
        // careful orientation of the ray coordinate ysstem above.
        let max_width = max(
            lerp(self.u_min, self.common.width[0], self.common.width[1]),
            lerp(self.u_max, self.common.width[0], self.common.width[1]),
        );
        if max(max(cp[0].y, cp[1].y), max(cp[2].y, cp[3].y)) + 0.5 * max_width < 0.0
            || min(min(cp[0].y, cp[1].y), min(cp[2].y, cp[3].y)) - 0.5 * max_width > 0.0
        {
            return None;
        }

        // Check for non-overlap in x.
        if max(max(cp[0].x, cp[1].x), max(cp[2].x, cp[3].x)) + 0.5 * max_width < 0.0
            || min(min(cp[0].x, cp[1].x), min(cp[2].x, cp[3].x)) - 0.5 * max_width > 0.0
        {
            return None;
        }

        // Check for non-overlap in z.
        let ray_length = ray.d.length();
        let z_max = ray_length * ray.t_max;
        if max(max(cp[0].z, cp[1].z), max(cp[2].z, cp[3].z)) + 0.5 * max_width < 0.0
            || min(min(cp[0].z, cp[1].z), min(cp[2].z, cp[3].z)) - 0.5 * max_width > z_max
        {
            return None;
        }

        // Compute refinement depth for curve, _max_depth_
        let l0 = (0..2).fold(0.0, |l, i| {
            max(
                l,
                max(
                    max(
                        abs(cp[i].x - 2.0 * cp[i + 1].x + cp[i + 2].x),
                        abs(cp[i].y - 2.0 * cp[i + 1].y + cp[i + 2].y),
                    ),
                    abs(cp[i].z - 2.0 * cp[i + 1].z + cp[i + 2].z),
                ),
            )
        });

        let eps = max(self.common.width[0], self.common.width[1]) * 0.05; // width / 20

        // Compute log base 4 by dividing log2 in half.
        let r0 = log2(1.41421356237 * 6.0 * l0 / (8.0 * eps)) / 2;
        let max_depth = clamp(r0, 0, 10);

        self.recursive_intersect(
            &ray,
            &cp,
            Arc::new(object_to_ray.inverse()),
            self.u_min,
            self.u_max,
            max_depth,
            None,
        )
    }

    /// Returns the surface area of the shape in object space.
    fn area(&self) -> Float {
        // Compute object-space control points for curve segment, cp_obj
        let cp_obj = self.blossom_bezier();
        let width0 = lerp(self.u_min, self.common.width[0], self.common.width[1]);
        let width1 = lerp(self.u_max, self.common.width[0], self.common.width[1]);
        let avg_width = (width0 + width1) * 0.5;
        let approx_length = (0..3).fold(0.0, |a, i| a + cp_obj[i].distance(cp_obj[i + 1]));
        approx_length * avg_width
    }
}

/// Common curve parameters
#[derive(Copy, Clone, Debug)]
pub struct CurveCommon {
    /// The curve type.
    pub curve_type: CurveType,

    /// Object space control points.
    pub cp_obj: [Point3f; 4],

    /// The width of the curve at the start and end points.
    pub width: [Float; 2],

    /// Surface normal at the start and end points.
    pub n: [Normal3f; 2],

    /// Angle between the two normal vectors.
    pub normal_angle: Float,

    /// 1 / sin(normal_angle).
    pub inv_sin_normal_angle: Float,
}

/// Create common parameters for curve.
///
/// The curve type.
/// * `curve_type` - Curve type.
/// * `c`          - Object space control points.
/// * `width`      - The width of the curve at the start and end points.
/// * `norm`       - Surface normal at the start and end points.
pub fn curve_common(
    curve_type: CurveType,
    c: [Point3f; 4],
    width: [Float; 2],
    norm: [Normal3f; 2],
) -> CurveCommon {
    let n = if norm.len() == 2 {
        [norm[0].normalize(), norm[1].normalize()]
    } else {
        [Normal3f::default(), Normal3f::default()]
    };

    let normal_angle = clamp(n[0].dot(&n[1]), 0.0, 1.0).acos();
    let inv_sin_normal_angle = 1.0 / normal_angle;

    CurveCommon {
        curve_type,
        cp_obj: c,
        n,
        width,
        normal_angle,
        inv_sin_normal_angle,
    }
}

/// Computes the blossom p(u0, u1, u2) of a cubic Bézier spline.
///
/// * `p`  - Control points.
/// * `u0` - The first u-extent.
/// * `u1` - The second u-extent.
/// * `u2` - The third u-extent.
fn blossom_bezier(p: &[Point3f], u0: Float, u1: Float, u2: Float) -> Point3f {
    let a = [
        lerp(u0, p[0], p[1]),
        lerp(u0, p[1], p[2]),
        lerp(u0, p[2], p[3]),
    ];
    let b = [lerp(u1, a[0], a[1]), lerp(u1, a[1], a[2])];
    lerp(u2, b[0], b[1])
}

/// Subdivides a Bézier curve and returns 7-control points; points 0 - 3 are
/// control points for first half of the split curve and points 3 - 6 for the
/// second half of the curve.
///
/// * `cp` - The control points.
fn subdivide_bezier(cp: &[Point3f]) -> [Point3f; 7] {
    [
        cp[0],
        (cp[0] + cp[1]) / 2.0,
        (cp[0] + 2.0 * cp[1] + cp[2]) / 4.0,
        (cp[0] + 3.0 * cp[1] + 3.0 * cp[2] + cp[3]) / 8.0,
        (cp[1] + 2.0 * cp[2] + cp[3]) / 4.0,
        (cp[2] + cp[3]) / 2.0,
        cp[3],
    ]
}

/// Evaluate a Bézier curve at given parameter and return the point and derivative
/// at the point.
///
/// * `cp` - The control points.
/// * `u`  - The parameter to evaluate.
fn eval_bezier(cp: &[Point3f], u: Float) -> (Point3f, Vector3f) {
    let cp1 = [
        lerp(u, cp[0], cp[1]),
        lerp(u, cp[1], cp[2]),
        lerp(u, cp[2], cp[3]),
    ];
    let cp2 = [lerp(u, cp1[0], cp1[1]), lerp(u, cp1[1], cp1[2])];

    let deriv = if (cp2[1] - cp2[0]).length_squared() > 0.0 {
        3.0 * (cp2[1] - cp2[0])
    } else {
        // For a cubic Bezier, if the first three control points (say) are
        // coincident, then the derivative of the curve is legitimately (0,0,0)
        // at u=0.  This is problematic for us, though, since we'd like to be
        // able to compute a surface normal there.  In that case, just punt and
        // take the difference between the first and last control points, which
        // ain't great, but will hopefully do.
        cp[3] - cp[0]
    };

    (lerp(u, cp2[0], cp2[1]), deriv)
}
