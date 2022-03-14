//! Curves

#![allow(dead_code)]
use core::geometry::*;
use core::interaction::*;
use core::paramset::*;
use core::pbrt::*;
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
    pub data: Arc<ShapeData>,

    /// Common curve parameters.
    pub common: CurveData,

    /// Minimum u-parameter for the curve.
    pub u_min: Float,

    /// Maximum u-parameter for the curve.
    pub u_max: Float,
}

impl Curve {
    /// Create a new curve.
    ///
    /// * `object_to_world`     - The object to world transfomation.
    /// * `world_to_object`     - The world to object transfomation.
    /// * `reverse_orientation` - Indicates whether their surface normal directions
    ///                           should be reversed from the default
    /// * `u_min`               - Minimum u-parameter for the curve.
    /// * `u_max`               - Maximum u-parameter for the curve.
    pub fn new(
        object_to_world: ArcTransform,
        world_to_object: ArcTransform,
        reverse_orientation: bool,
        common: CurveData,
        u_min: Float,
        u_max: Float,
    ) -> Self {
        Self {
            data: Arc::new(ShapeData::new(
                Arc::clone(&object_to_world),
                Some(Arc::clone(&world_to_object)),
                reverse_orientation,
            )),
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
    pub fn create_segments(
        object_to_world: ArcTransform,
        world_to_object: ArcTransform,
        reverse_orientation: bool,
        curve_type: CurveType,
        c: [Point3f; 4],
        width: [Float; 2],
        norm: Option<&[Normal3f]>,
        split_depth: i32, // Should really be usize.
    ) -> Vec<ArcShape> {
        let common = CurveData::new(curve_type, c, width, norm);

        let num_segments = 1 << split_depth;
        let mut segments = Vec::<ArcShape>::with_capacity(num_segments);

        let f = 1.0 / num_segments as Float;
        for i in 0..num_segments {
            let u_min = i as Float * f;
            let u_max = (i + 1) as Float * f;
            let curve = Curve::new(
                Arc::clone(&object_to_world),
                Arc::clone(&world_to_object),
                reverse_orientation,
                common,
                u_min,
                u_max,
            );
            segments.push(Arc::new(curve));
        }

        segments
    }

    /// Create `Curve`s from given parameter set, object to world transform,
    /// world to object transform and whether or not surface normal orientation
    /// is reversed.
    ///
    /// NOTE: Because we return a set of curves as `Vec<Arc<Shape>>` we cannot
    /// implement this as `From` trait :(
    ///
    /// * `p` - A tuple containing the parameter set, object to world transform,
    ///         world to object transform and whether or not surface normal
    ///         orientation is reversed.
    pub fn from_props(p: (&ParamSet, ArcTransform, ArcTransform, bool)) -> Vec<ArcShape> {
        let (params, o2w, w2o, reverse_orientation) = p;

        let width = params.find_one_float("width", 1.0);
        let width0 = params.find_one_float("width0", width);
        let width1 = params.find_one_float("width1", width);

        let degree = params.find_one_int("degree", 3_i32) as usize;
        if degree != 2 && degree != 3 {
            panic!(
                "Invalid degree {}: only degree 2 and 3 curves are supported.",
                degree
            );
        }

        let basis = params.find_one_string("basis", String::from("bezier"));
        if basis != "bezier" && basis != "bspline" {
            panic!(
                "Invalid basis '{}': only ''bezier' and 'bspline' are supported.",
                basis
            );
        }

        let cp = params.find_point3f("P");
        let ncp = cp.len();
        let n_segments: usize;
        if basis == "bezier" {
            // After the first segment, which uses degree+1 control points,
            // subsequent segments reuse the last control point of the previous
            // one and then use degree more control points.
            if ((ncp - 1 - degree) % degree) != 0 {
                panic!(
                    "Invalid number of control points {}: for the degree {} 
                    Bezier basis {} + n * {} are required, for n >= 0.",
                    ncp,
                    degree,
                    degree + 1,
                    degree
                );
            }
            n_segments = (ncp - 1) / degree;
        } else {
            if ncp < degree + 1 {
                panic!(
                    "Invalid number of control points {}: for the degree {} 
                      b-spline basis, must have >= {}.",
                    ncp,
                    degree,
                    degree + 1
                );
            }
            n_segments = ncp - degree;
        }

        let ctype = params.find_one_string("type", String::from("flat"));
        let curve_type = match &ctype[..] {
            "flat" => CurveType::Flat,
            "ribbon" => CurveType::Ribbon,
            "cylinder" => CurveType::Cylinder,
            t => {
                warn!("Unknown curve type '{}'.  Using 'cylinder'.", t);
                CurveType::Cylinder
            }
        };

        let mut n = params.find_normal3f("N");
        let nnorm = n.len();
        if nnorm > 0 {
            if curve_type != CurveType::Ribbon {
                warn!("Curve normals are only used with 'ribbon' type curves.");
                n = vec![];
            } else if nnorm != n_segments + 1 {
                panic!(
                    "Invalid number of normals {}: must provide {} normals for ribbon 
                    curves with {} segments.",
                    nnorm,
                    n_segments + 1,
                    n_segments
                );
            }
        } else if curve_type == CurveType::Ribbon {
            panic!("Must provide normals 'N' at curve endpoints with ribbon curves.");
        }

        let split_depth = params.find_one_float("splitdepth", 3.0) as i32;
        let sd = params.find_one_int("splitdepth", split_depth);

        let mut curves: Vec<ArcShape> = vec![];
        // Pointer to the first control point for the current segment. This is
        // updated after each loop iteration depending on the current basis.
        let mut cp_base = 0;
        for seg in 0..n_segments {
            let mut seg_cp_bezier = [Point3f::ZERO; 4];

            // First, compute the cubic Bezier control points for the current
            // segment and store them in segCpBezier. (It is admittedly
            // wasteful storage-wise to turn b-splines into Bezier segments and
            // wasteful computationally to turn quadratic curves into cubics,
            // but yolo.)
            if basis == "bezier" {
                if degree == 2 {
                    // Elevate to degree 3.
                    seg_cp_bezier[0] = cp[cp_base + 0];
                    seg_cp_bezier[1] = lerp(2.0 / 3.0, cp[cp_base + 0], cp[cp_base + 1]);
                    seg_cp_bezier[2] = lerp(1.0 / 3.0, cp[cp_base + 1], cp[cp_base + 2]);
                    seg_cp_bezier[3] = cp[cp_base + 2];
                } else {
                    // Allset.
                    for i in 0..4 {
                        seg_cp_bezier[i] = cp[cp_base + i];
                    }
                }
                cp_base += degree;
            } else {
                // Uniform b-spline.
                if degree == 2 {
                    // First compute equivalent Bezier control points via some
                    // blossiming.  We have three control points and a uniform
                    // knot vector; we'll label the points p01, p12, and p23.
                    // We want the Bezier control points of the equivalent
                    // curve, which are p11, p12, and p22.
                    let p01 = cp[cp_base + 0];
                    let p12 = cp[cp_base + 1];
                    let p23 = cp[cp_base + 2];

                    // We already have p12.
                    let p11 = lerp(0.5, p01, p12);
                    let p22 = lerp(0.5, p12, p23);

                    // Now elevate to degree 3.
                    seg_cp_bezier[0] = p11;
                    seg_cp_bezier[1] = lerp(2.0 / 3.0, p11, p12);
                    seg_cp_bezier[2] = lerp(1.0 / 3.0, p12, p22);
                    seg_cp_bezier[3] = p22;
                } else {
                    // Otherwise we will blossom from p011, p123, p234, and p345
                    // to the Bezier control points p222, p223, p233, and p333.
                    // https://people.eecs.berkeley.edu/~sequin/CS284/IMGS/cubicbsplinepoints.gif
                    let p012 = cp[cp_base + 0];
                    let p123 = cp[cp_base + 1];
                    let p234 = cp[cp_base + 2];
                    let p345 = cp[cp_base + 3];

                    let p122 = lerp(2.0 / 3.0, p012, p123);
                    let p223 = lerp(1.0 / 3.0, p123, p234);
                    let p233 = lerp(2.0 / 3.0, p123, p234);
                    let p334 = lerp(1.0 / 3.0, p234, p345);

                    let p222 = lerp(0.5, p122, p223);
                    let p333 = lerp(0.5, p233, p334);

                    seg_cp_bezier[0] = p222;
                    seg_cp_bezier[1] = p223;
                    seg_cp_bezier[2] = p233;
                    seg_cp_bezier[3] = p333;
                }
                cp_base += 1;
            }

            let width = [
                lerp(seg as Float / n_segments as Float, width0, width1),
                lerp((seg + 1) as Float / n_segments as Float, width0, width1),
            ];
            let c = Curve::create_segments(
                Arc::clone(&o2w),
                Arc::clone(&w2o),
                reverse_orientation,
                curve_type,
                seg_cp_bezier,
                width,
                if n.len() > 0 {
                    Some(&n[seg..seg + 2])
                } else {
                    None
                },
                sd,
            );
            curves.extend(c);
        }
        curves
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
    /// * `is_shadow_ray` - Used to terminate recursion on first hit for shadow rays.
    fn recursive_intersect<'scene, 'arena>(
        &self,
        ray: &Ray,
        cp: &[Point3f; 4],
        ray_to_object: &Transform,
        u0: Float,
        u1: Float,
        depth: u32,
        is_shadow_ray: bool,
    ) -> Option<Intersection<'scene>> {
        let ray_length = ray.d.length();

        if depth > 0 {
            // Split curve segment into sub-segments and test for intersection.
            let cp_split = subdivide_bezier(cp);

            // For each of the two segments, see if the ray's bounding box
            // overlaps the segment before recursively checking for
            // intersection with it.
            let u = [u0, (u0 + u1) / 2.0, u1];
            let mut hit: Option<Intersection<'scene>> = None;
            let mut cps = 0;
            for seg in 0..2 {
                // Splice containing the 4 control poitns for the current segment.
                //let cps = &cp_split[seg * 3..seg * 3 + 4];

                let max_width = max(
                    lerp(u[seg], self.common.width[0], self.common.width[1]),
                    lerp(u[seg + 1], self.common.width[0], self.common.width[1]),
                );

                // As above, check y first, since it most commonly lets us exit
                // out early.
                if max(
                    max(cp_split[cps + 0].y, cp_split[cps + 1].y),
                    max(cp_split[cps + 2].y, cp_split[cps + 3].y),
                ) + 0.5 * max_width
                    < 0.0
                    || min(
                        min(cp_split[cps + 0].y, cp_split[cps + 1].y),
                        min(cp_split[cps + 2].y, cp_split[cps + 3].y),
                    ) - 0.5 * max_width
                        > 0.0
                {
                    continue;
                }

                if max(
                    max(cp_split[cps + 0].x, cp_split[cps + 1].x),
                    max(cp_split[cps + 2].x, cp_split[cps + 3].x),
                ) + 0.5 * max_width
                    < 0.0
                    || min(
                        min(cp_split[cps + 0].x, cp_split[cps + 1].x),
                        min(cp_split[cps + 2].x, cp_split[cps + 3].x),
                    ) - 0.5 * max_width
                        > 0.0
                {
                    continue;
                }

                let z_max = ray_length * ray.t_max;
                if max(
                    max(cp_split[cps + 0].z, cp_split[cps + 1].z),
                    max(cp_split[cps + 2].z, cp_split[cps + 3].z),
                ) + 0.5 * max_width
                    < 0.0
                    || min(
                        min(cp_split[cps + 0].z, cp_split[cps + 1].z),
                        min(cp_split[cps + 2].z, cp_split[cps + 3].z),
                    ) - 0.5 * max_width
                        > z_max
                {
                    continue;
                }

                hit = self.recursive_intersect(
                    ray,
                    &[
                        cp_split[cps + 0],
                        cp_split[cps + 1],
                        cp_split[cps + 2],
                        cp_split[cps + 3],
                    ],
                    ray_to_object,
                    u[seg],
                    u[seg + 1],
                    depth - 1,
                    is_shadow_ray,
                );

                // If we found an intersection and this is a shadow ray,
                // we can exit out immediately.
                if hit.is_some() && is_shadow_ray {
                    return hit;
                }

                cps += 3;
            }
            hit
        } else {
            // Intersect ray with curve segment.

            // Test ray against segment endpoint boundaries.

            // Test sample point against tangent perpendicular at curve start.
            let mut edge = (cp[1].y - cp[0].y) * -cp[0].y + cp[0].x * (cp[0].x - cp[1].x);
            if edge < 0.0 {
                return None;
            }

            // Test sample point against tangent perpendicular at curve end.
            edge = (cp[2].y - cp[3].y) * -cp[3].y + cp[3].x * (cp[3].x - cp[2].x);
            if edge < 0.0 {
                return None;
            }

            // Compute line w that gives minimum distance to sample point.
            let segment_direction = Point2f::from(cp[3]) - Point2f::from(cp[0]);
            let denom = segment_direction.length_squared();
            if denom == 0.0 {
                return None;
            }
            let w = (-Vector2f::from(cp[0])).dot(&segment_direction) / denom;

            // Compute u coordinate of curve intersection point and hit_width.
            let u = clamp(lerp(w, u0, u1), u0, u1);
            let mut hit_width = lerp(u, self.common.width[0], self.common.width[1]);
            let mut n_hit = Normal3f::ZERO;
            if self.common.curve_type == CurveType::Ribbon {
                // Scale hit_width based on ribbon orientation.
                let sin0 =
                    sin((1.0 - u) * self.common.normal_angle) * self.common.inv_sin_normal_angle;
                let sin1 = sin(u * self.common.normal_angle) * self.common.inv_sin_normal_angle;
                n_hit = sin0 * self.common.n[0] + sin1 * self.common.n[1];
                hit_width *= n_hit.abs_dot(&ray.d) / ray_length;
            }

            // Test intersection point against curve width.
            let (pc, dpcdw) = eval_bezier(cp, clamp(w, 0.0, 1.0));
            let pt_curve_dist2 = pc.x * pc.x + pc.y * pc.y;
            if pt_curve_dist2 > hit_width * hit_width * 0.25 {
                return None;
            }
            let z_max = ray_length * ray.t_max;
            if pc.z < 0.0 || pc.z > z_max {
                return None;
            }

            // Compute v coordinate of curve intersection point.
            let pt_curve_dist = pt_curve_dist2.sqrt();
            let edge_func = dpcdw.x * (-pc.y) + pc.x * dpcdw.y;
            let v = if edge_func > 0.0 {
                0.5 + pt_curve_dist / hit_width
            } else {
                0.5 - pt_curve_dist / hit_width
            };

            // Compute hit `t` and partial derivatives for curve intersection.
            let t_hit = pc.z / ray_length;

            // Compute error bounds for curve intersection.
            let p_error = Vector3::new(2.0 * hit_width, 2.0 * hit_width, 2.0 * hit_width);

            // Compute dpdu and dpdv for curve intersection.
            let (_, dpdu) = eval_bezier(&self.common.cp_obj, u);
            assert!(
                dpdu != Vector3f::ZERO,
                "u={}, cp=[{:?}]",
                u,
                self.common.cp_obj
            );

            let dpdv = if self.common.curve_type == CurveType::Ribbon {
                Vector3::from(n_hit).cross(&dpdu).normalize() * hit_width
            } else {
                // Compute curve dpdv for flat and cylinder curves.
                let dpdu_plane = ray_to_object.inverse().transform_vector(&dpdu);
                let mut dpdv_plane =
                    Vector3::new(-dpdu_plane.y, dpdu_plane.x, 0.0).normalize() * hit_width;
                if self.common.curve_type == CurveType::Cylinder {
                    // Rotate dpdv_plane to give cylindrical appearance.
                    let theta = lerp(v, -90.0, 90.0);
                    let rot = Transform::rotate_axis(-theta, &dpdu_plane);
                    dpdv_plane = rot.transform_vector(&dpdv_plane);
                }
                ray_to_object.transform_vector(&dpdv_plane)
            };

            let mut si = SurfaceInteraction::new(
                ray.at(t_hit),
                p_error,
                Point2f::new(u, v),
                -ray.d,
                dpdu,
                dpdv,
                Normal3f::ZERO,
                Normal3f::ZERO,
                ray.time,
                Arc::clone(&self.data),
                0,
            );
            self.data
                .object_to_world
                .transform_surface_interaction(&mut si);

            Some(Intersection::new(t_hit, si))
        }
    }

    /// Returns geometric details if a ray intersects the shape intersection.
    /// If there is no intersection, `None` is returned.
    ///
    /// * `r`             - The ray.
    /// * `is_shadow_ray` - Used to terminate recursion on first hit for shadow rays.
    fn intersect<'scene, 'arena>(
        &self,
        r: &Ray,
        is_shadow_ray: bool,
    ) -> Option<Intersection<'scene>> {
        // Transform ray to object space.
        //
        // We could just use transform_ray() but there is minor adjustment in
        // it that adjusts t_max which is not in transform_ray_with_error().
        let (ray, _o_err, _d_err) = self
            .data
            .world_to_object
            .as_ref()
            .map(|w2o| w2o.transform_ray_with_error(r))
            .unwrap();

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
            let (dx_new, _dy) = coordinate_system(&ray.d);
            dx = dx_new;
        }

        let object_to_ray = Transform::look_at(&ray.o, &(ray.o + ray.d), &dx);
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

        // Compute refinement depth for curve, `max_depth`.
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

        // Compute log base 4 by dividing log2int in half.
        let r0 = (1.41421356237 * 6.0 * l0 / (8.0 * eps)).log2int() / 2;
        let max_depth = clamp(r0, 0, 10);

        self.recursive_intersect(
            &ray,
            &cp,
            &object_to_ray.inverse(),
            self.u_min,
            self.u_max,
            max_depth,
            is_shadow_ray,
        )
    }
}

impl Shape for Curve {
    /// Returns the shape type. Usually these are behind ArcShape and harder to
    /// debug. So this will be helpful.
    fn get_type(&self) -> &'static str {
        "curve"
    }

    /// Returns the underlying shape data.
    fn get_data(&self) -> Arc<ShapeData> {
        Arc::clone(&self.data)
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

        Bounds3f::new(cp_obj[0], cp_obj[1])
            .union(&Bounds3f::new(cp_obj[2], cp_obj[3]))
            .expand(max(width[0], width[1]) * 0.5)
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
        Self::intersect(self, r, false)
    }

    /// Returns `true` if a ray-shape intersection succeeds; otherwise `false`.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests; default to true.
    fn intersect_p(&self, r: &Ray, _test_alpha_texture: bool) -> bool {
        Self::intersect(self, r, true).is_some()
    }

    /// Returns the surface area of the shape in object space.
    fn area(&self) -> Float {
        // Compute object-space control points for curve segment, cp_obj.
        let cp_obj = self.blossom_bezier();
        let width0 = lerp(self.u_min, self.common.width[0], self.common.width[1]);
        let width1 = lerp(self.u_max, self.common.width[0], self.common.width[1]);
        let avg_width = (width0 + width1) * 0.5;
        let approx_length = (0..3).fold(0.0, |a, i| a + cp_obj[i].distance(cp_obj[i + 1]));
        approx_length * avg_width
    }

    /// Sample a point on the surface and return the PDF with respect to area on
    /// the surface.
    ///
    /// NOTE: The returned `Hit` value will have `wo` = Vector3f::ZERO.
    ///
    /// * `u` - Sample value to use.
    fn sample(&self, _u: &Point2f) -> (Hit, Float) {
        todo!()
    }
}

/// Common curve parameters
#[derive(Copy, Clone, Debug)]
pub struct CurveData {
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

impl CurveData {
    /// Create common parameters for curve.
    ///
    /// The curve type.
    /// * `curve_type` - Curve type.
    /// * `c`          - Object space control points.
    /// * `width`      - The width of the curve at the start and end points.
    /// * `norm`       - Surface normal at the start and end points.
    pub fn new(
        curve_type: CurveType,
        c: [Point3f; 4],
        width: [Float; 2],
        norm: Option<&[Normal3f]>,
    ) -> Self {
        let n = match norm {
            Some([n1, n2]) => [n1.normalize(), n2.normalize()],
            _ => [Normal3f::ZERO, Normal3f::ZERO],
        };

        let normal_angle = clamp(n[0].dot(&n[1]), 0.0, 1.0).acos();
        let inv_sin_normal_angle = 1.0 / normal_angle.sin();

        Self {
            curve_type,
            cp_obj: c,
            n,
            width,
            normal_angle,
            inv_sin_normal_angle,
        }
    }
}

/// Computes the blossom p(u0, u1, u2) of a cubic Bézier spline.
///
/// * `p`  - Control points.
/// * `u0` - The first u-extent.
/// * `u1` - The second u-extent.
/// * `u2` - The third u-extent.
fn blossom_bezier(p: &[Point3f; 4], u0: Float, u1: Float, u2: Float) -> Point3f {
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
fn subdivide_bezier(cp: &[Point3f; 4]) -> [Point3f; 7] {
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
fn eval_bezier(cp: &[Point3f; 4], u: Float) -> (Point3f, Vector3f) {
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
        // at u=0. This is problematic for us, though, since we'd like to be
        // able to compute a surface normal there. In that case, just punt and
        // take the difference between the first and last control points, which
        // ain't great, but will hopefully do.
        cp[3] - cp[0]
    };

    (lerp(u, cp2[0], cp2[1]), deriv)
}
