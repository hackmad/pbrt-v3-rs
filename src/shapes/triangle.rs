//! Triangles and triangle meshes

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;
use crate::core::sampling::*;
use crate::core::texture::*;
use crate::textures::*;
use std::collections::HashMap;
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

    /// Optional shadow alpha mask texture.
    pub shadow_alpha_mask: Option<ArcTexture<Float>>,

    /// Face indices.
    pub face_indices: Vec<usize>,
}

impl TriangleMesh {
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
    /// * `shadow_alpha_mask`   - Optional shadow alpha mask texture.
    /// * `face_indices`        - Face indices.
    pub fn new(
        object_to_world: ArcTransform,
        reverse_orientation: bool,
        vertex_indices: Vec<usize>,
        p: Vec<Point3f>,
        n: Vec<Normal3f>,
        s: Vec<Vector3f>,
        uv: Vec<Point2f>,
        alpha_mask: Option<ArcTexture<Float>>,
        shadow_alpha_mask: Option<ArcTexture<Float>>,
        face_indices: Vec<usize>,
    ) -> Self {
        let num_triangles = vertex_indices.len() % 3;
        assert!(num_triangles == 0);

        // Transform mesh vertices to world space.
        let tp = p.iter().map(|v| object_to_world.transform_point(&v));

        // Transform normals to world space.
        let tn = n.iter().map(|v| object_to_world.transform_normal(&v));

        // Transform normals to world space.
        let ts = s.iter().map(|v| object_to_world.transform_vector(&v));

        Self {
            num_triangles,
            vertex_indices,
            p: tp.collect(),
            n: tn.collect(),
            s: ts.collect(),
            uv,
            alpha_mask,
            shadow_alpha_mask,
            face_indices,
            data: ShapeData::new(object_to_world.clone(), None, reverse_orientation),
        }
    }

    /// Create a triangle mesh from vertex positions, normals, tangents, uv-coordinates
    /// and alpha mask.
    ///
    /// Returns a list of triangle data referencing it. Useful for shapes that
    /// convert to triangles.
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
    /// * `shadow_alpha_mask`   - Optional shadow alpha mask texture.
    /// * `face_indices`        - Face indices.
    pub fn create(
        object_to_world: ArcTransform,
        world_to_object: ArcTransform,
        reverse_orientation: bool,
        vertex_indices: Vec<usize>,
        p: Vec<Point3f>,
        n: Vec<Normal3f>,
        s: Vec<Vector3f>,
        uv: Vec<Point2f>,
        alpha_mask: Option<ArcTexture<Float>>,
        shadow_alpha_mask: Option<ArcTexture<Float>>,
        face_indices: Vec<usize>,
    ) -> Vec<ArcShape> {
        let num_triangles = vertex_indices.len() % 3;
        assert!(num_triangles == 0);

        let mesh = Self::new(
            object_to_world.clone(),
            reverse_orientation,
            vertex_indices,
            p,
            n,
            s,
            uv,
            alpha_mask,
            shadow_alpha_mask,
            face_indices,
        );

        let m = Arc::new(mesh);
        let mut tris = Vec::<ArcShape>::with_capacity(num_triangles);
        for i in 0..num_triangles {
            let tri = Triangle::new(
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

    /// Create a triangel mesh from given parameter set, object to world transform,
    /// world to object transform and whether or not surface normal orientation
    /// is reversed.
    ///
    /// NOTE: Because we return a set of curves as `Vec<Arc<Shape>>` we cannot
    /// implement this as `From` trait :(
    ///
    /// * `p`              - A tuple containing the parameter set, object to
    ///                      world transform, world to object transform and
    ///                      whether or not surface normal orientation is reversed.
    /// * `float_textures` - Float textures.
    pub fn from_props(
        p: (&ParamSet, ArcTransform, ArcTransform, bool),
        float_textures: &HashMap<String, ArcTexture<Float>>,
    ) -> Vec<ArcShape> {
        let (params, o2w, w2o, reverse_orientation) = p;

        let vi: Vec<usize> = params
            .find_int("indices")
            .iter()
            .map(|i| *i as usize)
            .collect();
        let nvi = vi.len();

        let p = params.find_point3f("P");
        let npi = p.len();

        let mut uvs = params.find_point2f("uv");
        if uvs.len() == 0 {
            uvs = params.find_point2f("st");
        }
        let mut nuvi = uvs.len();

        let mut temp_uvs: Vec<Point2f> = vec![];
        if uvs.len() == 0 {
            let mut fuv = params.find_float("uv");
            if fuv.len() == 0 {
                fuv = params.find_float("st");
            }
            nuvi = fuv.len();
            if nuvi > 0 {
                nuvi /= 2;
                for i in 0..nuvi {
                    temp_uvs.push(Point2f::new(fuv[2 * i], fuv[2 * i + 1]));
                }
                uvs = temp_uvs;
            }
        }
        if nuvi > 0 {
            if nuvi < npi {
                error!(
                    "Not enough of 'uv' for triangle mesh.  Expected {}, 
                    found {}.  Discarding.",
                    npi, nuvi
                );
                uvs = vec![];
            } else if nuvi > npi {
                error!(
                    "More 'uv' provided than will be used for triangle 
                    mesh.  ({} expcted, {} found)",
                    npi, nuvi
                );
            }
        }

        if nvi == 0 {
            error!("Vertex indices 'indices' not provided with triangle mesh shape");
            return vec![];
        }
        if npi == 0 {
            error!("Vertex positions 'P' not provided with triangle mesh shape");
            return vec![];
        }

        let mut s = params.find_vector3f("S");
        let nsi = s.len();
        if nsi > 0 && nsi != npi {
            error!("Number of 'S' for triangle mesh must match 'P'.");
            s = vec![];
        }

        let mut n = params.find_normal3f("N");
        let nni = n.len();
        if nni > 0 && nni != npi {
            error!("Number of 'N' for triangle mesh must match 'P'.");
            n = vec![];
        }
        for i in 0..nvi {
            if vi[i] >= npi {
                error!(
                    "trianglemesh has out-of-bounds vertex index {} ({} 'P' 
                    values were given",
                    vi[i], npi
                );
                return vec![];
            }
        }

        let mut face_indices: Vec<usize> = params
            .find_int("faceIndices")
            .iter()
            .map(|i| *i as usize)
            .collect();
        let nfi = face_indices.len();
        if nfi > 0 && nfi != nvi / 3 {
            error!(
                "Number of face indices, {}, doesn't match number of faces, {}",
                nfi,
                nvi / 3
            );
            face_indices = vec![];
        }

        let alpha_tex_name = params.find_one_texture("alpha", String::from(""));
        let alpha_tex = if alpha_tex_name.len() > 0 {
            if let Some(tex) = float_textures.get(&alpha_tex_name) {
                tex.clone()
            } else {
                warn!(
                    "Couldn't find float texture '{}' for 'alpha' parameter. 
                    Using float 'alpha' parameterer instead.",
                    alpha_tex_name
                );
                let alpha = params.find_one_float("alpha", 1.0);
                Arc::new(ConstantTexture::new(alpha))
            }
        } else {
            let alpha = params.find_one_float("alpha", 1.0);
            Arc::new(ConstantTexture::new(alpha))
        };

        let shadow_alpha_tex_name = params.find_one_texture("shadowalpha", String::from(""));
        let shadow_alpha_tex = if shadow_alpha_tex_name.len() > 0 {
            if let Some(tex) = float_textures.get(&shadow_alpha_tex_name) {
                tex.clone()
            } else {
                warn!(
                    "Couldn't find float texture '{}' for 'shadowalpha' 
                    parameter.  Using float 'shadowalpha' parameterer instead.",
                    alpha_tex_name
                );
                let alpha = params.find_one_float("shadowalpha", 1.0);
                Arc::new(ConstantTexture::new(alpha))
            }
        } else {
            let alpha = params.find_one_float("shadowalpha", 1.0);
            Arc::new(ConstantTexture::new(alpha))
        };

        Self::create(
            o2w.clone(),
            w2o.clone(),
            reverse_orientation,
            vi,
            p,
            n,
            s,
            uvs,
            Some(alpha_tex),
            Some(shadow_alpha_tex),
            face_indices,
        )
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

impl Triangle {
    /// Create a new triangle.
    ///
    /// * `object_to_world`     - The object to world transfomation.
    /// * `world_to_object`     - The world to object transfomation.
    /// * `reverse_orientation` - Indicates whether their surface normal directions
    ///                           should be reversed from the default
    /// * `mesh`                - The triangle mesh.
    /// * `triangle_index`      - The index of the triangle.
    pub fn new(
        object_to_world: ArcTransform,
        world_to_object: ArcTransform,
        reverse_orientation: bool,
        mesh: Arc<TriangleMesh>,
        triangle_index: usize,
    ) -> Self {
        Self {
            mesh: mesh.clone(),
            v: 3 * triangle_index,
            data: ShapeData::new(
                object_to_world.clone(),
                Some(world_to_object.clone()),
                reverse_orientation,
            ),
        }
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
            [
                Point2f::new(0.0, 0.0),
                Point2f::new(1.0, 0.0),
                Point2f::new(1.0, 1.0),
            ]
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
        (0..3).fold(Bounds3f::empty(), |b, i| {
            b.union(&world_to_object.transform_point(&self.mesh.p[self.v + i]))
        })
    }

    /// Returns a bounding box in the world space.
    ///
    /// Default is to transform the object bounds with the object-to0world
    /// transformation. Override for tighter bounds implementation.
    fn world_bound(&self) -> Bounds3f {
        (0..3).fold(Bounds3f::empty(), |b, i| b.union(&self.mesh.p[self.v + i]))
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
        let max_z_t = Vector3::new(p0t.z, p1t.z, p2t.z).abs().max_component();
        let delta_z = gamma(3) * max_z_t;

        // Compute delta_x and delta_y terms for triangle `t` error bounds.
        let max_x_t = Vector3::new(p0t.x, p1t.x, p2t.x).abs().max_component();
        let max_y_t = Vector3::new(p0t.y, p1t.y, p2t.y).abs().max_component();
        let delta_x = gamma(5) * (max_x_t + max_z_t);
        let delta_y = gamma(5) * (max_y_t + max_z_t);

        // Compute delta_e term for triangle `t` error bounds
        let delta_e = 2.0 * (gamma(2) * max_x_t * max_y_t + delta_y * max_x_t + delta_x * max_y_t);

        // Compute delta_t term for triangle `t` error bounds and check `t`.
        let max_e = Vector3::new(e0, e1, e2).abs().max_component();
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
            let (dpdu_new, dpdv_new) = coordinate_system(&ng.normalize());
            dpdu = dpdu_new;
            dpdv = dpdv_new;
        }

        // Compute error bounds for triangle intersection.
        let x_abs_sum = (b0 * p0.x).abs() + (b1 * p1.x).abs() + (b2 * p2.x).abs();
        let y_abs_sum = (b0 * p0.y).abs() + (b1 * p1.y).abs() + (b2 * p2.y).abs();
        let z_abs_sum = (b0 * p0.z).abs() + (b1 * p1.z).abs() + (b2 * p2.z).abs();
        let p_error = gamma(7) * Vector3::new(x_abs_sum, y_abs_sum, z_abs_sum);

        // Interpolate (u,v) parametric coordinates and hit point.
        let p_hit = b0 * p0 + b1 * p1 + b2 * p2;
        let uv_hit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

        // Test intersection against alpha texture, if present.
        if test_alpha_texture && !self.mesh.alpha_mask.is_none() {
            let isect_local = SurfaceInteraction::new(
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
        let mut isect = SurfaceInteraction::new(
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
                let (ss_new, ts_new) = coordinate_system(&ns.into());
                ss = ss_new;
                ts = ts_new;
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
                        let (dndu2, dndv2) = coordinate_system(&dn);
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

            isect.set_shading_geometry(ss, ts, dndu, dndv, true);
        }

        Some(Intersection::new(t, isect))
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
        let max_z_t = Vector3::new(p0t.z, p1t.z, p2t.z).abs().max_component();
        let delta_z = gamma(3) * max_z_t;

        // Compute delta_x and delta_y terms for triangle `t` error bounds.
        let max_x_t = Vector3::new(p0t.x, p1t.x, p2t.x).abs().max_component();
        let max_y_t = Vector3::new(p0t.y, p1t.y, p2t.y).abs().max_component();
        let delta_x = gamma(5) * (max_x_t + max_z_t);
        let delta_y = gamma(5) * (max_y_t + max_z_t);

        // Compute delta_e term for triangle `t` error bounds
        let delta_e = 2.0 * (gamma(2) * max_x_t * max_y_t + delta_y * max_x_t + delta_x * max_y_t);

        // Compute delta_t term for triangle `t` error bounds and check `t`.
        let max_e = Vector3::new(e0, e1, e2).abs().max_component();
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
                let (dpdu_new, dpdv_new) = coordinate_system(&ng.normalize());
                dpdu = dpdu_new;
                dpdv = dpdv_new;
            }

            // Interpolate parametric coordinates and hit point.
            let p_hit = b0 * p0 + b1 * p1 + b2 * p2;
            let uv_hit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

            let isect_local = SurfaceInteraction::new(
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

    /// Sample a point on the surface and return the PDF with respect to area on
    /// the surface.
    ///
    /// NOTE: The returned `Hit` value will have `wo` = Vector3f::default().
    ///
    /// * `u` - Sample value to use.
    fn sample_area(&self, u: &Point2f) -> (Hit, Float) {
        let b = uniform_sample_triangle(u);

        // Get triangle vertices in `p0`, `p1`, and `p2`.
        let p0 = self.mesh.p[self.v];
        let p1 = self.mesh.p[self.v + 1];
        let p2 = self.mesh.p[self.v + 2];

        let p = b[0] * p0 + b[1] * p1 + (1.0 - b[0] - b[1]) * p2;

        // Compute surface normal for sampled point on triangle.
        let mut n = Normal3f::from((p1 - p0).cross(&(p2 - p0))).normalize();

        // Ensure correct orientation of the geometric normal; follow the same
        // approach as was used in intersect().
        if self.mesh.n.len() > 0 {
            let ns = Vector3f::from(
                b[0] * self.mesh.n[self.v]
                    + b[1] * self.mesh.n[self.v + 1]
                    + (1.0 - b[0] - b[1]) * self.mesh.n[self.v + 2],
            );
            n = n.face_forward(&ns);
        } else if self.data.reverse_orientation ^ self.data.transform_swaps_handedness {
            n *= -1.0;
        }

        // Compute error bounds for sampled point on triangle.
        let p_abs_sum = (b[0] * p0).abs() + (b[1] * p1).abs() + ((1.0 - b[0] - b[1]) * p2).abs();
        let p_error = gamma(6) * Vector3f::new(p_abs_sum.x, p_abs_sum.y, p_abs_sum.z);
        let it = Hit::new(p, 0.0, p_error, Vector3f::default(), n, None);
        let pdf = 1.0 / self.area();
        (it, pdf)
    }

    /// Sample a point on the shape given a reference point and return the PDF
    /// with respect to the solid angle from ref.
    ///
    /// * `hit` - Reference point on shape.
    /// * `u`   - Sample value to use.
    fn sample_solid_angle(&self, _hit: &Hit, _u: &Point2f) -> (Hit, Float) {
        todo!()
    }

    /// Returns the PDF with respect to solid angle.
    ///
    /// * `hit` - The interaction hit point.
    /// * `wi`  - The incoming direction.
    fn pdf_solid_angle(&self, _hit: &Hit, _wi: &Vector3f) -> Float {
        todo!()
    }
}
