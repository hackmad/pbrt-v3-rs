//! Triangles and triangle meshes

#![allow(dead_code)]
use super::{
    coordinate_system, empty_bounds3, gamma, intersection, point2, shape_data, surface_interaction,
    vector3, ArcShape, ArcTexture, ArcTransform, Bounds3f, Float, Intersection, Normal3, Normal3f,
    Point2f, Point3f, Ray, Shape, ShapeData, Union, Vector3, Vector3f,
};
use std::mem::size_of;
use std::sync::Arc;

/// Triangle mesh
#[derive(Clone)]
pub struct TriangleMesh {
    /// Common shape data.
    pub data: ShapeData,

    /// The number of triangles.
    pub num_triangles: usize,

    /// Vertex indices. For the ith triangle, its three vertex positions are
    /// p[vertex_indices[3 * i]], p[vertex_indices[3 * i + 1]], and
    /// p[vertex_indices[3 * i + 2]]
    pub vertex_indices: Vec<usize>,

    /// Vertex positions.
    pub p: Vec<Point3f>,

    /// Vertex normals. This will be empty if there are none.
    pub n: Vec<Normal3f>,

    /// Tangent vectors per vertex. This will be empty if there are none.
    pub s: Vec<Vector3f>,

    /// Paramteric uv-coordinates per vertex. This will be empty if there are none.
    pub uv: Vec<Point2f>,

    /// Optional alpha mask texture, which can be used to cut away parts of
    /// triangle surfaces
    pub alpha_mask: Option<ArcTexture<Float>>,
}

/// Create a new triangle mesh.
///
/// * `object_to_world`     - The object to world transfomation.
/// * `reverse_orientation` - Indicates whether their surface normal directions
///                           should be reversed from the default
/// * `vertex_indices`      - Vertex indices for triangles. For the ith triangle,
///                           its three vertex positions are p[vertex_indices[3 * i]],
/// * `p`                   - Vertex positions.
///                           p[vertex_indices[3 * i + 1]], and
///                           p[vertex_indices[3 * i + 2]]
/// * `n`                   - Vertex normals.
/// * `s`                   - Tangent vectors per vertex.
/// * `uv`                  - Paramteric uv-coordinates.
/// * `alpha_mask`          - Optional alpha mask texture, which can be used to
///                           cut away parts of triangle surfaces
pub fn triangle_mesh(
    object_to_world: ArcTransform,
    reverse_orientation: bool,
    vertex_indices: Vec<usize>,
    p: Vec<Point3f>,
    n: Vec<Normal3f>,
    s: Vec<Vector3f>,
    uv: Vec<Point2f>,
    alpha_mask: Option<ArcTexture<Float>>,
) -> TriangleMesh {
    let num_triangles = vertex_indices.len() % 3;
    assert!(num_triangles == 0);

    // Transform mesh vertices to world space.
    let tp = p.iter().map(|v| object_to_world.transform_point(&v));

    // Transform normals to world space.
    let tn = n.iter().map(|v| object_to_world.transform_normal(&v));

    // Transform normals to world space.
    let ts = s.iter().map(|v| object_to_world.transform_vector(&v));

    TriangleMesh {
        num_triangles,
        vertex_indices,
        p: tp.collect(),
        n: tn.collect(),
        s: ts.collect(),
        uv,
        alpha_mask,
        data: shape_data(object_to_world.clone(), None, reverse_orientation),
    }
}

/// Triangle.
#[derive(Clone)]
pub struct Triangle {
    /// Common shape data.
    pub data: ShapeData,

    /// The mesh.
    pub mesh: Arc<TriangleMesh>,

    /// The index of the first vertex in the triangle. The other two are
    /// calculated as v + 1 and v + 2.
    pub v: usize,
}

/// Create a new triangle.
///
/// * `object_to_world`     - The object to world transfomation.
/// * `world_to_object`     - The world to object transfomation.
/// * `reverse_orientation` - Indicates whether their surface normal directions
///                           should be reversed from the default
/// * `mesh`                - The triangle mesh.
/// * `triangle_index`      - The index of the triangle.
pub fn triangle(
    object_to_world: ArcTransform,
    world_to_object: ArcTransform,
    reverse_orientation: bool,
    mesh: Arc<TriangleMesh>,
    triangle_index: usize,
) -> Triangle {
    Triangle {
        mesh: mesh.clone(),
        v: 3 * triangle_index,
        data: shape_data(
            object_to_world.clone(),
            Some(world_to_object.clone()),
            reverse_orientation,
        ),
    }
}

impl Triangle {
    /// Returns the uv-coordinates for the triangle. If there are no uv
    /// coordinates, then default ones [(0,0), (1,0), (1,1)] are returned.
    fn get_uvs(&self) -> [Point2f; 3] {
        if self.mesh.uv.len() > 0 {
            [
                self.mesh.uv[self.v],
                self.mesh.uv[self.v + 1],
                self.mesh.uv[self.v + 2],
            ]
        } else {
            [point2(0.0, 0.0), point2(1.0, 0.0), point2(1.0, 1.0)]
        }
    }
}

impl Shape for Triangle {
    /// Returns the underlying shape data.
    fn get_data(&self) -> ShapeData {
        self.data.clone()
    }

    /// Returns a bounding box in the shapes object space.
    fn object_bound(&self) -> Bounds3f {
        // We can unwrap safely because the factory methods guarantee world_to_object
        // is passed. If it is constructed without that, then tough luck!
        let world_to_object = self.data.world_to_object.clone().unwrap();
        (0..3).fold(empty_bounds3::<Float>(), |b, i| {
            b.union(&world_to_object.transform_point(&self.mesh.p[self.v + i]))
        })
    }

    /// Returns a bounding box in the world space.
    ///
    /// Default is to transform the object bounds with the object-to0world
    /// transformation. Override for tighter bounds implementation.
    fn world_bound(&self) -> Bounds3f {
        (0..3).fold(empty_bounds3::<Float>(), |b, i| {
            b.union(&self.mesh.p[self.v + i])
        })
    }

    /// Returns geometric details if a ray intersects the shape intersection.
    /// If there is no intersection, `None` is returned.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests.
    fn intersect<'a>(&self, r: &Ray, test_alpha_texture: bool) -> Option<Intersection<'a>> {
        // Get triangle vertices in p0, p1, and p2
        let p0 = self.mesh.p[self.v];
        let p1 = self.mesh.p[self.v + 1];
        let p2 = self.mesh.p[self.v + 2];

        // Perform ray-triangle intersection test.

        // Transform triangle vertices based on ray origin.

        // Translate vertices based on ray origin.
        let o = Vector3::from(r.o);
        let mut p0t = p0 - o;
        let mut p1t = p1 - o;
        let mut p2t = p2 - o;

        // Permute components of triangle vertices and ray direction.
        let kz = r.d.abs().max_dimension();
        let kx = kz + 1;
        let ky = kx + 1;
        let d = r.d.permute(kx, ky, kz);
        p0t = p0t.permute(kx, ky, kz);
        p1t = p1t.permute(kx, ky, kz);
        p2t = p2t.permute(kx, ky, kz);

        // Apply shear transformation to translated vertex positions.
        let sx = -d.x / d.z;
        let sy = -d.y / d.z;
        let sz = 1.0 / d.z;
        p0t.x += sx * p0t.z;
        p0t.y += sy * p0t.z;
        p1t.x += sx * p1t.z;
        p1t.y += sy * p1t.z;
        p2t.x += sx * p2t.z;
        p2t.y += sy * p2t.z;

        // Compute edge function coefficients e0, e1, e2.
        let mut e0 = p1t.x * p2t.y - p1t.y * p2t.x;
        let mut e1 = p2t.x * p0t.y - p2t.y * p0t.x;
        let mut e2 = p0t.x * p1t.y - p0t.y * p1t.x;

        // Fallback to double-precision test at triangle edges.
        if size_of::<Float>() == size_of::<f32>() && (e0 == 0.0 || e1 == 0.0 || e2 == 0.0) {
            let p2txp1ty = (p2t.x as f64) * (p1t.y as f64);
            let p2typ1tx = (p2t.y as f64) * (p1t.x as f64);
            e0 = (p2typ1tx - p2txp1ty) as Float;

            let p0txp2ty = (p0t.x as f64) * (p2t.y as f64);
            let p0typ2tx = (p0t.y as f64) * (p2t.x as f64);
            e1 = (p0typ2tx - p0txp2ty) as Float;

            let p1txp0ty = (p1t.x as f64) * (p0t.y as f64);
            let p1typ0tx = (p1t.y as f64) * (p0t.x as f64);
            e2 = (p1typ0tx - p1txp0ty) as Float;
        }

        // Perform triangle edge and determinant tests.
        if (e0 < 0.0 || e1 < 0.0 || e2 < 0.0) && (e0 > 0.0 || e1 > 0.0 || e2 > 0.0) {
            return None;
        }

        let det = e0 + e1 + e2;
        if det == 0.0 {
            return None;
        }

        // Compute scaled hit distance to triange and test against ray `t` range.
        p0t.z *= sz;
        p1t.z *= sz;
        p2t.z *= sz;
        let t_scaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        if det < 0.0 && (t_scaled >= 0.0 || t_scaled < r.t_max * det) {
            return None;
        } else if det > 0.0 && (t_scaled <= 0.0 || t_scaled > r.t_max * det) {
            return None;
        }

        // Compute barycentric coordinates and `t` value for triangle intersection.
        let inv_det = 1.0 / det;
        let b0 = e0 * inv_det;
        let b1 = e1 * inv_det;
        let b2 = e2 * inv_det;
        let t = t_scaled * inv_det;

        // Ensure that computed triangle is conservatively greater than zero.

        // Compute delta_z term for triangle `t` error bounds.
        let max_z_t = vector3(p0t.z, p1t.z, p2t.z).abs().max_component();
        let delta_z = gamma(3) * max_z_t;

        // Compute delta_x and delta_y terms for triangle `t` error bounds.
        let max_x_t = vector3(p0t.x, p1t.x, p2t.x).abs().max_component();
        let max_y_t = vector3(p0t.y, p1t.y, p2t.y).abs().max_component();
        let delta_x = gamma(5) * (max_x_t + max_z_t);
        let delta_y = gamma(5) * (max_y_t + max_z_t);

        // Compute delta_e term for triangle `t` error bounds
        let delta_e = 2.0 * (gamma(2) * max_x_t * max_y_t + delta_y * max_x_t + delta_x * max_y_t);

        // Compute delta_t term for triangle `t` error bounds and check `t`.
        let max_e = vector3(e0, e1, e2).abs().max_component();
        let delta_t = 3.0
            * (gamma(3) * max_e * max_z_t + delta_e * max_z_t + delta_z * max_e)
            * inv_det.abs();
        if t <= delta_t {
            return None;
        }

        // Compute triangle partial derivatives.
        let uv = self.get_uvs();

        // Compute deltas for triangle partial derivatives.
        let duv02 = uv[0] - uv[2];
        let duv12 = uv[1] - uv[2];
        let dp02 = p0 - p2;
        let dp12 = p1 - p2;
        let determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
        let degenerate_uv = determinant.abs() < 1e-8;
        let mut dpdu = Vector3f::default();
        let mut dpdv = Vector3f::default();
        if !degenerate_uv {
            let invdet = 1.0 / determinant;
            dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
            dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
        }
        if degenerate_uv || dpdu.cross(&dpdv).length_squared() == 0.0 {
            // Handle zero determinant for triangle partial derivative matrix.
            let ng = (p2 - p0).cross(&(p1 - p0));
            if ng.length_squared() == 0.0 {
                // The triangle is actually degenerate; the intersection is bogus.
                return None;
            }
            coordinate_system(&ng.normalize(), &mut dpdu, &mut dpdv);
        }

        // Compute error bounds for triangle intersection.
        let x_abs_sum = (b0 * p0.x).abs() + (b1 * p1.x).abs() + (b2 * p2.x).abs();
        let y_abs_sum = (b0 * p0.y).abs() + (b1 * p1.y).abs() + (b2 * p2.y).abs();
        let z_abs_sum = (b0 * p0.z).abs() + (b1 * p1.z).abs() + (b2 * p2.z).abs();
        let p_error = gamma(7) * vector3(x_abs_sum, y_abs_sum, z_abs_sum);

        // Interpolate (u,v) parametric coordinates and hit point.
        let p_hit = b0 * p0 + b1 * p1 + b2 * p2;
        let uv_hit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

        // Test intersection against alpha texture, if present.
        if test_alpha_texture && !self.mesh.alpha_mask.is_none() {
            let isect_local = surface_interaction(
                p_hit,
                Vector3f::default(),
                uv_hit,
                -r.d,
                dpdu,
                dpdv,
                Normal3f::default(),
                Normal3f::default(),
                r.time,
                Some(Arc::new(self.clone())),
            );

            let alpha_mask = self.mesh.alpha_mask.clone().unwrap();
            if alpha_mask.evaluate(&isect_local) == 0.0 {
                return None;
            }
        }

        // Fill SurfaceInteraction from triangle hit.
        let mut isect = surface_interaction(
            p_hit,
            p_error,
            uv_hit,
            -r.d,
            dpdu,
            dpdv,
            Normal3f::default(),
            Normal3f::default(),
            r.time,
            Some(Arc::new(self.clone())),
        );

        // Override surface normal in isect for triangle.
        isect.hit.n = Normal3::from(dp02.cross(&dp12).normalize());
        if self.get_data().reverse_orientation ^ self.get_data().transform_swaps_handedness {
            isect.hit.n = -isect.hit.n;
        }
        isect.shading.n = isect.hit.n;

        let has_vertex_normals = self.mesh.n.len() > 0;
        let has_vertex_tangents = self.mesh.s.len() > 0;
        if has_vertex_normals || has_vertex_tangents {
            // Initialize triangle shading geometry.

            // Compute shading normal ns for triangle.
            let mut ns = isect.hit.n;
            if has_vertex_normals {
                let ns2 = b0 * self.mesh.n[self.v]
                    + b1 * self.mesh.n[self.v + 1]
                    + b2 * self.mesh.n[self.v + 2];
                if ns2.length_squared() > 0.0 {
                    ns = ns2.normalize();
                }
            };

            // Compute shading tangent ss for triangle.
            let mut ss = isect.dpdu;
            if has_vertex_tangents {
                let ss2 = b0 * self.mesh.s[self.v]
                    + b1 * self.mesh.s[self.v + 1]
                    + b2 * self.mesh.s[self.v + 2];
                if ss2.length_squared() > 0.0 {
                    ss = ss2;
                }
            };
            ss = ss.normalize();

            // Compute shading bitangent ts for triangle and adjust ss
            let mut ts = ss.cross(&ns.into());
            if ts.length_squared() > 0.0 {
                ts = ts.normalize();
                ss = ts.cross(&ns.into());
            } else {
                coordinate_system(&ns.into(), &mut ss, &mut ts);
            }

            // Compute dndu and dndv for triangle shading geometry.
            let (dndu, dndv) = if has_vertex_normals {
                // Compute deltas for triangle partial derivatives of normal
                let duv02 = uv[0] - uv[2];
                let duv12 = uv[1] - uv[2];
                let dn1 = self.mesh.n[self.v] - self.mesh.n[self.v + 2];
                let dn2 = self.mesh.n[self.v + 1] - self.mesh.n[self.v + 2];

                let determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
                let degenerate_uv = determinant.abs() < 1e-8;
                if degenerate_uv {
                    // We can still compute dndu and dndv, with respect to the
                    // same arbitrary coordinate system we use to compute dpdu
                    // and dpdv when this happens. It's important to do this
                    // (rather than giving up) so that ray differentials for
                    // rays reflected from triangles with degenerate
                    // parameterizations are still reasonable.
                    let dn = Vector3::from(self.mesh.n[self.v + 2] - self.mesh.n[self.v]).cross(
                        &Vector3::from(self.mesh.n[self.v + 1] - self.mesh.n[self.v]),
                    );
                    if dn.length_squared() == 0.0 {
                        (Normal3f::default(), Normal3f::default())
                    } else {
                        let mut dndu2 = Vector3f::default();
                        let mut dndv2 = Vector3f::default();
                        coordinate_system(&dn, &mut dndu2, &mut dndv2);
                        (Normal3::from(dndu2), Normal3::from(dndv2))
                    }
                } else {
                    let invdet = 1.0 / determinant;
                    (
                        (duv12[1] * dn1 - duv02[1] * dn2) * invdet,
                        (-duv12[0] * dn1 + duv02[0] * dn2) * invdet,
                    )
                }
            } else {
                (Normal3f::default(), Normal3f::default())
            };

            if self.get_data().reverse_orientation {
                ts = -ts;
            }

            let (new_n, new_shading) = isect.update_shading_geometry(ss, ts, dndu, dndv, true);
            isect.hit.n = new_n;
            isect.shading = new_shading;
        }

        Some(intersection(t, isect))
    }

    /// Returns `true` if a ray-shape intersection succeeds; otherwise `false`.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests.
    fn intersect_p(&self, r: &Ray, test_alpha_texture: bool) -> bool {
        // Get triangle vertices in p0, p1, and p2
        let p0 = self.mesh.p[self.v];
        let p1 = self.mesh.p[self.v + 1];
        let p2 = self.mesh.p[self.v + 2];

        // Perform ray-triangle intersection test.

        // Transform triangle vertices based on ray origin.

        // Translate vertices based on ray origin.
        let o = Vector3::from(r.o);
        let mut p0t = p0 - o;
        let mut p1t = p1 - o;
        let mut p2t = p2 - o;

        // Permute components of triangle vertices and ray direction.
        let kz = r.d.abs().max_dimension();
        let kx = kz + 1;
        let ky = kx + 1;
        let d = r.d.permute(kx, ky, kz);
        p0t = p0t.permute(kx, ky, kz);
        p1t = p1t.permute(kx, ky, kz);
        p2t = p2t.permute(kx, ky, kz);

        // Apply shear transformation to translated vertex positions.
        let sx = -d.x / d.z;
        let sy = -d.y / d.z;
        let sz = 1.0 / d.z;
        p0t.x += sx * p0t.z;
        p0t.y += sy * p0t.z;
        p1t.x += sx * p1t.z;
        p1t.y += sy * p1t.z;
        p2t.x += sx * p2t.z;
        p2t.y += sy * p2t.z;

        // Compute edge function coefficients e0, e1, e2.
        let mut e0 = p1t.x * p2t.y - p1t.y * p2t.x;
        let mut e1 = p2t.x * p0t.y - p2t.y * p0t.x;
        let mut e2 = p0t.x * p1t.y - p0t.y * p1t.x;

        // Fallback to double-precision test at triangle edges.
        if size_of::<Float>() == size_of::<f32>() && (e0 == 0.0 || e1 == 0.0 || e2 == 0.0) {
            let p2txp1ty = (p2t.x as f64) * (p1t.y as f64);
            let p2typ1tx = (p2t.y as f64) * (p1t.x as f64);
            e0 = (p2typ1tx - p2txp1ty) as Float;

            let p0txp2ty = (p0t.x as f64) * (p2t.y as f64);
            let p0typ2tx = (p0t.y as f64) * (p2t.x as f64);
            e1 = (p0typ2tx - p0txp2ty) as Float;

            let p1txp0ty = (p1t.x as f64) * (p0t.y as f64);
            let p1typ0tx = (p1t.y as f64) * (p0t.x as f64);
            e2 = (p1typ0tx - p1txp0ty) as Float;
        }

        // Perform triangle edge and determinant tests.
        if (e0 < 0.0 || e1 < 0.0 || e2 < 0.0) && (e0 > 0.0 || e1 > 0.0 || e2 > 0.0) {
            return false;
        }

        let det = e0 + e1 + e2;
        if det == 0.0 {
            return false;
        }

        // Compute scaled hit distance to triange and test against ray `t` range.
        p0t.z *= sz;
        p1t.z *= sz;
        p2t.z *= sz;
        let t_scaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        if det < 0.0 && (t_scaled >= 0.0 || t_scaled < r.t_max * det) {
            return false;
        } else if det > 0.0 && (t_scaled <= 0.0 || t_scaled > r.t_max * det) {
            return false;
        }

        // Compute barycentric coordinates and value for triangle intersection.
        let inv_det = 1.0 / det;
        let b0 = e0 * inv_det;
        let b1 = e1 * inv_det;
        let b2 = e2 * inv_det;
        let t = t_scaled * inv_det;

        // Ensure that computed triangle is conservatively greater than zero.

        // Compute delta_z term for triangle `t` error bounds.
        let max_z_t = vector3(p0t.z, p1t.z, p2t.z).abs().max_component();
        let delta_z = gamma(3) * max_z_t;

        // Compute delta_x and delta_y terms for triangle `t` error bounds.
        let max_x_t = vector3(p0t.x, p1t.x, p2t.x).abs().max_component();
        let max_y_t = vector3(p0t.y, p1t.y, p2t.y).abs().max_component();
        let delta_x = gamma(5) * (max_x_t + max_z_t);
        let delta_y = gamma(5) * (max_y_t + max_z_t);

        // Compute delta_e term for triangle `t` error bounds
        let delta_e = 2.0 * (gamma(2) * max_x_t * max_y_t + delta_y * max_x_t + delta_x * max_y_t);

        // Compute delta_t term for triangle `t` error bounds and check `t`.
        let max_e = vector3(e0, e1, e2).abs().max_component();
        let delta_t = 3.0
            * (gamma(3) * max_e * max_z_t + delta_e * max_z_t + delta_z * max_e)
            * inv_det.abs();
        if t <= delta_t {
            return false;
        }

        // Test intersection against alpha texture, if present.
        if test_alpha_texture && !self.mesh.alpha_mask.is_none() {
            // Compute triangle partial derivatives.
            let uv = self.get_uvs();

            // Compute deltas for triangle partial derivatives.
            let duv02 = uv[0] - uv[2];
            let duv12 = uv[1] - uv[2];
            let dp02 = p0 - p2;
            let dp12 = p1 - p2;
            let determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
            let degenerate_uv = determinant.abs() < 1e-8;

            let mut dpdu = Vector3f::default();
            let mut dpdv = Vector3f::default();
            if !degenerate_uv {
                let invdet = 1.0 / determinant;
                dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
                dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
            }
            if degenerate_uv || dpdu.cross(&dpdv).length_squared() == 0.0 {
                // Handle zero determinant for triangle partial derivative matrix.
                let ng = (p2 - p0).cross(&(p1 - p0));
                if ng.length_squared() == 0.0 {
                    // The triangle is actually degenerate; the intersection is bogus.
                    return false;
                }
                coordinate_system(&ng.normalize(), &mut dpdu, &mut dpdv);
            }

            // Interpolate parametric coordinates and hit point.
            let p_hit = b0 * p0 + b1 * p1 + b2 * p2;
            let uv_hit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

            let isect_local = surface_interaction(
                p_hit,
                Vector3f::default(),
                uv_hit,
                -r.d,
                dpdu,
                dpdv,
                Normal3f::default(),
                Normal3f::default(),
                r.time,
                Some(Arc::new(self.clone())),
            );

            let alpha_mask = self.mesh.alpha_mask.clone().unwrap();
            if alpha_mask.evaluate(&isect_local) == 0.0 {
                return false;
            }
        }

        true
    }

    /// Returns the surface area of the shape in object space.
    fn area(&self) -> Float {
        let p0 = self.mesh.p[self.v];
        let p1 = self.mesh.p[self.v + 1];
        let p2 = self.mesh.p[self.v + 2];
        0.5 * (p1 - p0).cross(&(p2 - p0)).length()
    }
}

/// Create a triangle mesh and the return a list of triangle data referencing it.
/// Useful for shapes that convert to triangles.
///
/// * `object_to_world`     - The object to world transfomation.
/// * `reverse_orientation` - Indicates whether their surface normal directions
///                           should be reversed from the default
/// * `vertex_indices`      - Vertex indices for triangles. For the ith triangle,
///                           its three vertex positions are p[vertex_indices[3 * i]],
///                           p[vertex_indices[3 * i + 1]], and
///                           p[vertex_indices[3 * i + 2]]
/// * `p`                   - Vertex positions.
/// * `n`                   - Vertex normals.
/// * `s`                   - Tangent vectors per vertex.
/// * `uv`                  - Paramteric uv-coordinates.
/// * `alpha_mask`          - Optional alpha mask texture, which can be used to
///                           cut away parts of triangle surfaces
pub fn create_triangle_mesh(
    object_to_world: ArcTransform,
    world_to_object: ArcTransform,
    reverse_orientation: bool,
    vertex_indices: Vec<usize>,
    p: Vec<Point3f>,
    n: Vec<Normal3f>,
    s: Vec<Vector3f>,
    uv: Vec<Point2f>,
    alpha_mask: Option<ArcTexture<Float>>,
) -> Vec<ArcShape> {
    let num_triangles = vertex_indices.len() % 3;
    assert!(num_triangles == 0);

    let mesh = triangle_mesh(
        object_to_world.clone(),
        reverse_orientation,
        vertex_indices,
        p,
        n,
        s,
        uv,
        alpha_mask,
    );

    let m = Arc::new(mesh);
    let mut tris = Vec::<ArcShape>::with_capacity(num_triangles);
    for i in 0..num_triangles {
        let tri = triangle(
            object_to_world.clone(),
            world_to_object.clone(),
            reverse_orientation,
            m.clone(),
            i,
        );
        tris.push(Arc::new(tri));
    }

    tris
}