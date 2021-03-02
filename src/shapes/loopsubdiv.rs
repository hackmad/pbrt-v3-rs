//! Loop Subdivision Surfaces.

#![allow(dead_code)]

use super::TriangleMesh;
use crate::core::geometry::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;
use std::cmp::{Ordering, PartialOrd};
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::sync::Arc;

/// Subdivision surface vertex.
///
/// NOTE: We use i64 and -1 to indicate no value instead of using Option<usize>
/// to make the algorithm a little bit cleaner and clearer to implement following
/// the original PBRT implementation.
#[derive(Clone)]
struct SDVertex {
    /// Index in vertices.
    pub vi: i64,

    /// Vertex position.
    pub p: Point3f,

    /// An arbitrary face adjacent to this vertex.
    pub start_face: i64,

    /// First vertex in next level of subdivision.
    pub child: i64,

    /// Indicates regular or extraordinary vertex. Interior vertices with valence
    /// other than six, or boundary vertices with valence other than four, are
    /// called extraordinary vertices; otherwise, they are called regular.
    pub regular: bool,

    /// Indicates vertex is on the boundary of the mesh.
    pub boundary: bool,
}

impl SDVertex {
    /// Create a new subdivision vertex.
    ///
    /// * `vi` - Vertex index.
    /// * `p`  - Vertex position.
    pub fn new(vi: i64, p: Point3f) -> Self {
        Self {
            vi,
            p,
            start_face: -1,
            child: -1,
            regular: false,
            boundary: false,
        }
    }

    /// Returns number of vertices directly adjacent to this vertex.
    ///
    /// * `face` - The list of faces.
    pub fn valence(&self, faces: &[Arc<SDFace>]) -> usize {
        let mut f = self.start_face;
        let mut nf = 1_usize;
        if !self.boundary {
            while f != -1 {
                f = faces[f as usize].next_face(self.vi);
                if f == self.start_face {
                    break;
                }
                nf += 1;
            }
            nf
        } else {
            while f != -1 {
                f = faces[f as usize].next_face(self.vi);
                if f != -1 {
                    nf += 1;
                }
            }
            f = self.start_face;
            while f != -1 {
                f = faces[f as usize].prev_face(self.vi);
                if f != -1 {
                    nf += 1;
                }
            }
            nf + 1
        }
    }

    /// Returns position of vertices around this vertex.
    ///
    /// * `verts` - The list of vertices.
    /// * `faces` - The of faces.
    pub fn one_ring(&self, verts: &[Arc<SDVertex>], faces: &[Arc<SDFace>]) -> Vec<Point3f> {
        let mut p: Vec<Point3f> = vec![];
        let mut f = self.start_face;
        if !self.boundary {
            loop {
                let nv = faces[f as usize].next_vert(self.vi) as usize;
                p.push(verts[nv].p);
                f = faces[f as usize].next_face(self.vi);
                if f == self.start_face {
                    break;
                }
            }
        } else {
            let mut f2: i64;
            loop {
                f2 = faces[f as usize].next_face(self.vi);
                if f2 == -1 {
                    break;
                } else {
                    f = f2;
                }
            }
            let nv = faces[f as usize].next_vert(self.vi) as usize;
            p.push(verts[nv].p);
            loop {
                let nv = faces[f as usize].prev_vert(self.vi) as usize;
                p.push(verts[nv].p);
                f = faces[f as usize].prev_face(self.vi);
                if f == -1 {
                    break;
                }
            }
        }
        p
    }
}

impl Default for SDVertex {
    /// Return a default value for `SDVertex`.
    fn default() -> Self {
        Self {
            vi: -1,
            p: Point3f::default(),
            start_face: -1,
            child: -1,
            regular: false,
            boundary: false,
        }
    }
}

impl PartialEq for SDVertex {
    /// Determines equality based on the position.
    fn eq(&self, other: &Self) -> bool {
        self.p == other.p
    }
}

// Subdivision surface edge.
///
/// NOTE: We use i64 and -1 to indicate no value instead of using Option<usize>
/// to make the algorithm a little bit cleaner and clearer to implement following
/// the original PBRT implementation.
#[derive(Clone, Eq)]
struct SDEdge {
    /// Vertex indices for the endpoints of the edge.
    pub v: [i64; 2],

    /// The face indices of adjacent faces.
    pub f: [i64; 2],

    /// The index of first face found that is adjacent to this edge.
    pub f0_edge_num: i64,
}

impl SDEdge {
    /// Create a new subdvision surface edge.
    ///
    /// * `v0` - First endpoint.
    /// * `v1` - Second endpoint.
    fn new(v0: i64, v1: i64) -> Self {
        Self {
            v: [v0, v1],
            f: [-1, -1],
            f0_edge_num: -1,
        }
    }
}

impl Default for SDEdge {
    /// Return a default value for `SDEdge`.
    fn default() -> Self {
        Self {
            v: [-1, -1],
            f: [-1, -1],
            f0_edge_num: -1,
        }
    }
}

impl PartialEq for SDEdge {
    /// Determine equality based on vertex indices.
    ///
    /// * `other` - The `SDEdge` to compare to.
    fn eq(&self, other: &Self) -> bool {
        self.v[0] == other.v[0] && self.v[1] == other.v[1]
    }
}

impl PartialOrd for SDEdge {
    /// Define an ordering for edges based on vertex indices.
    ///
    /// * `other` - The other edge.
    fn partial_cmp(&self, other: &SDEdge) -> Option<Ordering> {
        if self.v[0] == other.v[0] {
            Some(self.v[1].cmp(&other.v[1]))
        } else {
            Some(self.v[0].cmp(&other.v[0]))
        }
    }
}

impl Hash for SDEdge {
    // Feeds the vertex indices into the given `Hasher`.
    //
    // * `state` - Hasher state.
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.v[0].hash(state);
        self.v[1].hash(state);
    }
}

/// Subdivision surface face.
///
/// NOTE: We use i64 and -1 to indicate no value instead of using Option<usize>
/// to make the algorithm a little bit cleaner and clearer to implement following
/// the original PBRT implementation.
#[derive(Clone)]
struct SDFace {
    /// Index of vertices of the face.
    pub v: [i64; 3],

    /// Adjacent face indices.
    pub f: [i64; 3],

    /// Indices of faces at the next level of subdivision.
    pub children: [i64; 4],
}

impl SDFace {
    /// Find the index in v[] of a given vertex.
    ///
    /// * `vert` - The vertex index.
    pub fn vnum(&self, vert: i64) -> i64 {
        for i in 0..3 {
            if self.v[i] == (vert as i64) {
                return i as i64;
            }
        }
        -1
    }

    /// Returns index of next face adjacent to given vertex.
    ///
    /// * `vert` - The vertex index.
    pub fn next_face(&self, vert: i64) -> i64 {
        let i = self.vnum(vert);
        if i == -1 {
            panic!("next_face({}) not found", vert);
        }
        self.f[i as usize]
    }

    /// Returns index of previous face adjacent to given vertex.
    ///
    /// * `vert` - The vertex index.
    pub fn prev_face(&self, vert: i64) -> i64 {
        let i = self.vnum(vert);
        if i == -1 {
            panic!("prev_face({}) not found", vert);
        }
        self.f[prev(i) as usize]
    }

    /// Returns index of next vertex adjacent to given vertex.
    ///
    /// * `vert` - The vertex index.
    pub fn next_vert(&self, vert: i64) -> i64 {
        let i = self.vnum(vert);
        if i == -1 {
            panic!("next_vert({}) not found", vert);
        }
        self.v[next(i) as usize]
    }

    /// Returns index of previous vertex adjacent to given vertex.
    ///
    /// * `vert` - The vertex index.
    pub fn prev_vert(&self, vert: i64) -> i64 {
        let i = self.vnum(vert);
        if i == -1 {
            panic!("prev_vert({}) not found", vert);
        }
        self.v[prev(i) as usize]
    }

    /// Returns the vertex index of vertex opposite an edge.
    ///
    /// * `v0` - Vertex index of first endpoint of edge.
    /// * `v1` - Vertex index of second endpoint of edge.
    pub fn other_vert(&self, v0: i64, v1: i64) -> i64 {
        for i in 0..3 {
            if self.v[i] != v0 && self.v[i] != v1 {
                return self.v[i];
            }
        }
        panic!("Basic logic error in SDVertex::other_vert({}, {})", v0, v1);
    }
}

impl Default for SDFace {
    /// Return a default value for `SDFace`.
    fn default() -> Self {
        Self {
            v: [-1, -1, -1],
            f: [-1, -1, -1],
            children: [-1, -1, -1, -1],
        }
    }
}

/// Implements Loop Subdivision.
pub struct LoopSubDiv {}

impl LoopSubDiv {
    /// Subdivide a triangle mesh using loop subdivision.
    /// * `object_to_world`     - The object to world transfomation.
    /// * `world_to_object`     - The world to object transfomation.
    /// * `reverse_orientation` - Indicates whether their surface normal directions
    ///                           should be reversed from the default
    /// * `n_levels`            - Number of subdivision levels.
    /// * `vertex_indices`      - Vertex indices.
    /// * `p`                   - Vertex positions.
    pub fn subdivide(
        object_to_world: ArcTransform,
        world_to_object: ArcTransform,
        reverse_orientation: bool,
        n_levels: usize,
        vertex_indices: Vec<usize>,
        p: Vec<Point3f>,
    ) -> Vec<ArcShape> {
        // Allocate `LoopSubDiv` vertices and faces.
        let n_vertices = p.len();
        let mut vertices: Vec<Arc<SDVertex>> = Vec::with_capacity(n_vertices);

        for i in 0..n_vertices {
            vertices.push(Arc::new(SDVertex::new(i as i64, p[i])));
        }

        let n_faces = vertex_indices.len() / 3;
        let mut faces: Vec<Arc<SDFace>> = Vec::with_capacity(n_faces);
        for _i in 0..n_faces {
            faces.push(Arc::new(SDFace::default()));
        }

        // Set face to vertex pointers.
        let mut vp = &vertex_indices[0..];
        for i in 0..n_faces {
            for j in 0..3 {
                let v = Arc::get_mut(&mut vertices[vp[j]]).unwrap();
                v.start_face = i as i64;

                let f = Arc::get_mut(&mut faces[i]).unwrap();
                f.v[j] = v.vi;
            }
            vp = &vp[3..];
        }

        // Set neighbor pointers in `faces`.
        let mut edges: HashSet<SDEdge> = HashSet::new();
        for i in 0..n_faces {
            for edge_num in 0..3 {
                // Update neighbor pointer for `edge_num`.
                let (v0, v1) = (edge_num as usize, next(edge_num) as usize);
                let mut e = SDEdge::new(faces[i].v[v0], faces[i].v[v1]);
                if !edges.contains(&e) {
                    // Handle new edge.
                    e.f[0] = i as i64;
                    e.f0_edge_num = edge_num;
                    edges.insert(e);
                } else {
                    // Handle previously seen edge.
                    let e = edges.take(&e).unwrap();
                    let f = Arc::get_mut(&mut faces[i]).unwrap();
                    f.f[e.f0_edge_num as usize] = i as i64;
                    f.f[edge_num as usize] = e.f[0];
                }
            }
        }

        // Finish vertex initialization.
        for i in 0..n_vertices {
            let v = Arc::get_mut(&mut vertices[i]).unwrap();
            let mut f = v.start_face;
            loop {
                f = faces[f as usize].next_face(v.vi);
                if f == -1 || f == v.start_face {
                    break;
                }
            }
            let valence = v.valence(&faces);
            v.boundary = f == -1;
            if !v.boundary && valence == 6 {
                v.regular = true;
            } else if v.boundary && valence == 4 {
                v.regular = true;
            } else {
                v.regular = false;
            }
        }

        // Refine `LoopSubDiv` into triangles.
        let mut f = faces;
        let mut v = vertices;

        for _i in 0..n_levels {
            // Update `f` and `v` for next level of subdivision.
            let mut new_faces: Vec<Arc<SDFace>> = vec![];
            let mut new_vertices: Vec<Arc<SDVertex>> = vec![];

            // Allocate next level of children in mesh tree.
            for vertex in v.iter_mut() {
                let vi = new_vertices.len();
                Arc::get_mut(vertex).unwrap().child = vi as i64;

                let mut child = SDVertex::default();
                child.vi = vi as i64;
                child.regular = vertex.regular;
                child.boundary = vertex.boundary;
                new_vertices.push(Arc::new(child));
            }
            for face in f.iter_mut() {
                for k in 0..4 {
                    let fi = new_faces.len();
                    Arc::get_mut(face).unwrap().children[k] = fi as i64;
                    new_faces.push(Arc::new(SDFace::default()));
                }
            }

            // Update vertex positions and create new edge vertices.

            // Update vertex positions for even vertices.
            for vi in 0..v.len() {
                let ci = v[vi].child as usize;
                if let Some(child) = Arc::get_mut(&mut new_vertices[ci]) {
                    if !v[vi].boundary {
                        // Apply one-ring rule for even vertex.
                        if v[vi].regular {
                            child.p = weight_one_ring(v[vi].clone(), 1.0 / 16.0, &v, &f);
                        } else {
                            child.p =
                                weight_one_ring(v[vi].clone(), beta(v[vi].valence(&f)), &v, &f);
                        }
                    } else {
                        // Apply boundary rule for even vertex.
                        child.p = weight_boundary(v[vi].clone(), 1.0 / 8.0, &v, &f);
                    }
                }
            }

            // Compute new odd edge vertices.
            let mut edge_verts: HashMap<SDEdge, i64> = HashMap::new();
            for face in f.iter() {
                for k in 0..3 {
                    // Compute odd vertex on `k`th edge.
                    let edge = SDEdge::new(face.v[k], face.v[next(k as i64) as usize]);
                    if !edge_verts.contains_key(&edge) {
                        // Create and initialize new odd vertex
                        let mut vert = SDVertex::default();
                        vert.vi = new_vertices.len() as i64;
                        vert.regular = true;
                        vert.boundary = face.f[k] == -1;
                        vert.start_face = face.children[3];

                        // Apply edge rules to compute new vertex position.
                        let v0 = v[edge.v[0] as usize].clone();
                        let v1 = v[edge.v[1] as usize].clone();

                        if vert.boundary {
                            vert.p = 0.5 * v0.p + 0.5 * v1.p;
                        } else {
                            let ov = v[face.other_vert(v0.vi, v1.vi) as usize].clone();
                            let fk = f[face.f[k] as usize].clone();
                            let fk_ov = v[fk.other_vert(edge.v[0], edge.v[1]) as usize].clone();
                            vert.p = 3.0 / 8.0 * v0.p
                                + 3.0 / 8.0 * v1.p
                                + 1.0 / 8.0 * ov.p
                                + 1.0 / 8.0 * fk_ov.p;
                        }
                        edge_verts.insert(edge, vert.vi);

                        new_vertices.push(Arc::new(vert));
                    }
                }
            }

            // Update new mesh topology.

            // Update even vertex face pointers.
            for vertex in v.iter() {
                if vertex.child != -1 {
                    let start_face = f[vertex.start_face as usize].clone();
                    let vert_num = start_face.vnum(vertex.vi) as usize;
                    let child_start_face = start_face.children[vert_num];
                    if let Some(child) = Arc::get_mut(&mut new_vertices[vertex.child as usize]) {
                        child.start_face = child_start_face;
                    }
                }
            }

            // Update face neighbor pointers.
            for face in f.iter() {
                for j in 0..3_usize {
                    // Update children `f` pointers for siblings.
                    let c3 = face.children[3] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[c3]) {
                        child.f[j] = face.children[next(j as i64) as usize];
                    }
                    let cj = face.children[j] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[cj]) {
                        child.f[next(j as i64) as usize] = face.children[3];
                    }

                    // Update children `f` pointers for neighbor children.
                    let f2 = face.f[j];
                    if let Some(child) = Arc::get_mut(&mut new_faces[cj]) {
                        child.f[j] = if f2 != -1 {
                            let face2 = f[f2 as usize].clone();
                            face2.children[face2.vnum(face.v[j]) as usize]
                        } else {
                            -1
                        };
                    }

                    let f2 = face.f[prev(j as i64) as usize];
                    let cj = face.children[j] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[cj]) {
                        child.f[prev(j as i64) as usize] = if f2 != -1 {
                            let face2 = f[f2 as usize].clone();
                            face2.children[face2.vnum(face.v[j]) as usize]
                        } else {
                            -1
                        };
                    }
                }
            }

            // Update face vertex pointers
            for face in f.iter() {
                for j in 0..3_usize {
                    // Update child vertex pointer to new even vertex.
                    let cj = face.children[j] as usize;
                    let mut vi = -1_i64;
                    if let Some(child) = Arc::get_mut(&mut new_faces[cj]) {
                        child.v[j] = v[face.v[j] as usize].child;

                        // Update child vertex pointer to new odd vertex
                        let edge = SDEdge::new(face.v[j], face.v[next(j as i64) as usize]);
                        if let Some(vert) = edge_verts.get(&edge) {
                            vi = *vert;
                            child.v[next(j as i64) as usize] = vi;
                        }
                    }

                    let cnextj = face.children[next(j as i64) as usize] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[cnextj]) {
                        child.v[j] = vi;
                    }

                    let c3 = face.children[3] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[c3]) {
                        child.v[j] = vi;
                    }
                }
            }

            // Prepare for next level of subdivision.
            f = new_faces.split_off(0);
            v = new_vertices.split_off(0);
        }

        // Push vertices to limit surface.
        let mut p_limit: Vec<Point3f> = Vec::with_capacity(v.len());
        for i in 0..v.len() {
            if v[i].boundary {
                p_limit.push(weight_boundary(v[i].clone(), 1.0 / 5.0, &v, &f));
            } else {
                p_limit.push(weight_one_ring(
                    v[i].clone(),
                    loop_gamma(v[i].clone().valence(&f)),
                    &v,
                    &f,
                ));
            }
        }
        for i in 0..v.len() {
            Arc::get_mut(&mut v[i]).unwrap().p = p_limit[i];
        }

        // Compute vertex tangents on limit surface.
        let mut ns: Vec<Normal3f> = Vec::with_capacity(v.len());
        for i in 0..v.len() {
            let vertex = v[i].clone();
            let mut s = Vector3f::default();
            let mut t = Vector3f::default();
            let valence = vertex.valence(&f);
            let p_ring = vertex.one_ring(&v, &f);
            if !vertex.boundary {
                // Compute tangents of interior face.
                for j in 0..valence {
                    let theta = TWO_PI * j as Float / valence as Float;
                    s += cos(theta) * Vector3f::from(p_ring[j]);
                    t += sin(theta) * Vector3f::from(p_ring[j]);
                }
            } else {
                // Compute tangents of boundary face
                s = p_ring[valence - 1] - p_ring[0];
                if valence == 2 {
                    t = Vector3f::from(p_ring[0] + p_ring[1] - 2.0 * vertex.p);
                } else if valence == 3 {
                    t = p_ring[1] - vertex.p;
                } else if valence == 4 {
                    // regular
                    t = Vector3f::from(
                        -1.0 * p_ring[0]
                            + 2.0 * p_ring[1]
                            + 2.0 * p_ring[2]
                            + -1.0 * p_ring[3]
                            + -2.0 * vertex.p,
                    );
                } else {
                    let theta = PI / (valence - 1) as Float;
                    t = Vector3f::from(sin(theta) * (p_ring[0] + p_ring[valence - 1]));
                    for k in 1..valence - 1 {
                        let wt = (2.0 * cos(theta) - 2.0) * sin(k as Float * theta);
                        t += Vector3f::from(wt * p_ring[k]);
                    }
                    t = -t;
                }
            }
            ns.push(Normal3f::from(s.cross(&t)));
        }

        // Create triangle mesh from subdivision mesh.
        let n_tris = f.len();
        let mut vertex_indices: Vec<usize> = Vec::with_capacity(3 * n_tris);
        for face in f.iter() {
            for j in 0..3 {
                vertex_indices.push(face.v[j] as usize);
            }
        }

        TriangleMesh::create(
            object_to_world.clone(),
            world_to_object.clone(),
            reverse_orientation,
            vertex_indices,
            p_limit,
            ns,
            vec![],
            vec![],
            None,
            None,
            vec![],
        )
    }

    /// Create `LoopSubDiv` from given parameter set, object to world transform,
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

        let n_levels = params.find_one_int("nlevels", 3) as usize;
        let vertex_indices: Vec<usize> = params
            .find_int("indices")
            .iter()
            .map(|i| *i as usize)
            .collect();
        let p = params.find_point3f("P");
        if vertex_indices.len() == 0 {
            panic!("Vertex indices 'indices' not provided for LoopSubDiv shape.");
        }
        if p.len() == 0 {
            panic!("Vertex positions 'P' not provided for LoopSubDiv shape.");
        }

        Self::subdivide(
            o2w.clone(),
            w2o.clone(),
            reverse_orientation,
            n_levels,
            vertex_indices,
            p,
        )
    }
}

/// Compute the bea value based on the vertex's valence to ensure smoothness.
///
/// * `valence` - Valence of a vertex.
fn beta(valence: usize) -> Float {
    if valence == 3 {
        3.0 / 16.0
    } else {
        3.0 / (8.0 * valence as Float)
    }
}

/// Compute the appropriate vertex weights based on the valence of the vertex.
///
/// * `valence` - Valence of a vertex.
fn loop_gamma(valence: usize) -> Float {
    1.0 / (valence as Float + 3.0 / (8.0 * beta(valence)))
}

/// Applies given weight to the one-ring of adjacent vertices and returns a
/// new position.
///
/// * `vert`  - The vertex.
/// * `beta`  - The weight.
/// * `verts` - The list of vertices.
/// * `faces` - The list of faces.
fn weight_one_ring(
    vert: Arc<SDVertex>,
    beta: Float,
    verts: &[Arc<SDVertex>],
    faces: &[Arc<SDFace>],
) -> Point3f {
    let valence = vert.valence(faces);
    let p_ring = vert.one_ring(verts, faces);
    let mut p = vert.p * (1.0 as Float - valence as Float * beta);
    for i in 0..valence {
        p += p_ring[i] * beta;
    }
    p
}

/// Applies given weights at the boundary vertex and returns a new position.
///
/// * `vert`  - The vertex.
/// * `beta`  - The weight.
/// * `verts` - The list of vertices.
/// * `faces` - The list of faces.
fn weight_boundary(
    vert: Arc<SDVertex>,
    beta: Float,
    verts: &[Arc<SDVertex>],
    faces: &[Arc<SDFace>],
) -> Point3f {
    // put _vert_ one-ring in _p_ring_
    let valence = vert.valence(faces);
    let p_ring = vert.one_ring(verts, faces);
    let mut p = vert.p * (1.0 as Float - 2.0 as Float * beta);
    p += p_ring[0] * beta;
    p += p_ring[(valence - 1) as usize] * beta;
    p
}

/// Returns the next index modulo 3.
///
/// `i`  - The index.
fn next(i: i64) -> i64 {
    (i + 1) % 3
}

/// Returns the previous index modulo 3.
///
/// `i`  - The index.
fn prev(i: i64) -> i64 {
    (i + 2) % 3
}
