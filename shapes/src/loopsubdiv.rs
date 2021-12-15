//! Loop Subdivision Surfaces.

#![allow(dead_code)]

use super::TriangleMesh;
use core::geometry::*;
use core::paramset::*;
use core::pbrt::*;
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
    /// * `p`  - Vertex position.
    pub fn new(p: Point3f) -> Self {
        Self {
            p,
            start_face: -1,
            child: -1,
            regular: false,
            boundary: false,
        }
    }

    /// Returns number of vertices directly adjacent to this vertex.
    ///
    /// * `vi`   - Vertex index.
    /// * `face` - The list of faces.
    pub fn valence(&self, vi: i64, faces: &[Arc<SDFace>]) -> usize {
        let mut f = self.start_face;
        let mut nf = 1_usize;

        if !self.boundary {
            // Compute valence of interior vertex.
            while f != -1 {
                f = faces[f as usize].next_face(vi);
                if f == self.start_face {
                    break;
                }
                nf += 1;
            }
            nf
        } else {
            // Compute valence of boundary vertex.
            while f != -1 {
                f = faces[f as usize].next_face(vi);
                if f != -1 {
                    nf += 1;
                } else {
                    break;
                }
            }
            f = self.start_face;
            while f != -1 {
                f = faces[f as usize].prev_face(vi);
                if f != -1 {
                    nf += 1;
                } else {
                    break;
                }
            }
            nf + 1
        }
    }

    /// Returns position of vertices around this vertex.
    ///
    /// * `vi`    - The vertex index.
    /// * `verts` - The list of vertices.
    /// * `faces` - The of faces.
    pub fn one_ring(
        &self,
        vi: i64,
        verts: &[Arc<SDVertex>],
        faces: &[Arc<SDFace>],
    ) -> Vec<Point3f> {
        let mut p: Vec<Point3f> = vec![];
        if !self.boundary {
            let mut f = self.start_face;
            loop {
                let nv = faces[f as usize].next_vert(vi) as usize;
                p.push(verts[nv].p);
                f = faces[f as usize].next_face(vi);
                if f == self.start_face {
                    break;
                }
            }
        } else {
            let mut f = self.start_face;
            let mut f2: i64;
            loop {
                f2 = faces[f as usize].next_face(vi);
                if f2 == -1 {
                    break;
                } else {
                    f = f2;
                }
            }
            let nv = faces[f as usize].next_vert(vi) as usize;
            p.push(verts[nv].p);
            loop {
                let nv = faces[f as usize].prev_vert(vi) as usize;
                p.push(verts[nv].p);
                f = faces[f as usize].prev_face(vi);
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
            v: [min(v0, v1), max(v0, v1)],
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
        let mut verts: Vec<Arc<SDVertex>> =
            p.iter().map(|&pos| Arc::new(SDVertex::new(pos))).collect();

        let n_faces = vertex_indices.len() / 3;
        let mut faces: Vec<Arc<SDFace>> =
            (0..n_faces).map(|_| Arc::new(SDFace::default())).collect();

        // Set face to vertex pointers.
        for i in 0..n_faces {
            for j in 0..3 {
                let vi = vertex_indices[i * 3 + j];

                let f = Arc::get_mut(&mut faces[i]).unwrap();
                f.v[j] = vi as i64;

                let v = Arc::get_mut(&mut verts[vi]).unwrap();
                v.start_face = i as i64;
            }
        }

        // Set neighbor pointers in `faces`.
        let mut edges: HashSet<SDEdge> = HashSet::new();
        for i in 0..n_faces {
            for edge_num in 0..3 {
                // Update neighbour pointer for `edge_num`.
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
                    let f = Arc::get_mut(&mut faces[e.f[0] as usize]).unwrap();
                    f.f[e.f0_edge_num as usize] = i as i64;

                    let f = Arc::get_mut(&mut faces[i as usize]).unwrap();
                    f.f[edge_num as usize] = e.f[0];
                }
            }
        }

        // Finish vertex initialization.
        for i in 0..verts.len() {
            let v = Arc::get_mut(&mut verts[i]).unwrap();
            let mut f = v.start_face;
            loop {
                f = faces[f as usize].next_face(i as i64);
                if f == -1 || f == v.start_face {
                    break;
                }
            }
            let valence = v.valence(i as i64, &faces);
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
        for _i in 0..n_levels {
            // Update `faces` and `verts` for next level of subdivision.
            let mut new_faces: Vec<Arc<SDFace>> = vec![];
            let mut new_vertices: Vec<Arc<SDVertex>> = vec![];

            // Allocate next level of children in mesh tree.
            for vert in verts.iter_mut() {
                let vi = new_vertices.len();
                Arc::get_mut(vert).unwrap().child = vi as i64;

                let mut child = SDVertex::default();
                child.regular = vert.regular;
                child.boundary = vert.boundary;
                new_vertices.push(Arc::new(child));
            }
            for face in faces.iter_mut() {
                for k in 0..4 {
                    let fi = new_faces.len();
                    Arc::get_mut(face).unwrap().children[k] = fi as i64;
                    new_faces.push(Arc::new(SDFace::default()));
                }
            }

            // Update vertex positions and create new edge vertices.

            // Update vertex positions for even vertices.
            for vi in 0..verts.len() {
                let ci = verts[vi].child as usize;
                if let Some(child) = Arc::get_mut(&mut new_vertices[ci]) {
                    if !verts[vi].boundary {
                        // Apply one-ring rule for even vertex.
                        if verts[vi].regular {
                            child.p = weight_one_ring(
                                Arc::clone(&verts[vi]),
                                1.0 / 16.0,
                                vi as i64,
                                &verts,
                                &faces,
                            );
                        } else {
                            child.p = weight_one_ring(
                                Arc::clone(&verts[vi]),
                                beta(verts[vi].valence(vi as i64, &faces)),
                                vi as i64,
                                &verts,
                                &faces,
                            );
                        }
                    } else {
                        // Apply boundary rule for even vertex.
                        child.p = weight_boundary(
                            Arc::clone(&verts[vi]),
                            1.0 / 8.0,
                            vi as i64,
                            &verts,
                            &faces,
                        );
                    }
                }
            }

            // Compute new odd edge vertices.
            let mut edge_verts: HashMap<SDEdge, i64> = HashMap::new();
            for i in 0..faces.len() {
                for k in 0..3 {
                    // Compute odd vertex on `k`th edge.
                    let edge = SDEdge::new(faces[i].v[k], faces[i].v[next(k as i64) as usize]);
                    if !edge_verts.contains_key(&edge) {
                        // Create and initialize new odd vertex
                        let new_vi = new_vertices.len();
                        let mut vert = SDVertex::default();
                        vert.regular = true;
                        vert.boundary = faces[i].f[k] == -1;
                        vert.start_face = faces[i].children[3];

                        // Apply edge rules to compute new vertex position.
                        let v0 = Arc::clone(&verts[edge.v[0] as usize]);
                        let v1 = Arc::clone(&verts[edge.v[1] as usize]);

                        if vert.boundary {
                            vert.p = 0.5 * v0.p + 0.5 * v1.p;
                        } else {
                            let ov = Arc::clone(
                                &verts[faces[i].other_vert(edge.v[0], edge.v[1]) as usize],
                            );
                            let fk = Arc::clone(&faces[faces[i].f[k] as usize]);
                            let fk_ov =
                                Arc::clone(&verts[fk.other_vert(edge.v[0], edge.v[1]) as usize]);
                            vert.p = 3.0 / 8.0 * v0.p
                                + 3.0 / 8.0 * v1.p
                                + 1.0 / 8.0 * ov.p
                                + 1.0 / 8.0 * fk_ov.p;
                        }
                        edge_verts.insert(edge, new_vi as i64);

                        new_vertices.push(Arc::new(vert));
                    }
                }
            }

            // Update new mesh topology.

            // Update even vertex face pointers.
            for (vi, vert) in verts.iter().enumerate() {
                if vert.child != -1 {
                    let start_face = Arc::clone(&faces[vert.start_face as usize]);
                    let vert_num = start_face.vnum(vi as i64) as usize;
                    let child_start_face = start_face.children[vert_num];
                    if let Some(child) = Arc::get_mut(&mut new_vertices[vert.child as usize]) {
                        child.start_face = child_start_face;
                    }
                }
            }

            // Update face neighbor pointers.
            for i in 0..faces.len() {
                for j in 0..3_usize {
                    // Update children `f` pointers for siblings.
                    let c3 = faces[i].children[3] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[c3]) {
                        child.f[j] = faces[i].children[next(j as i64) as usize];
                    }
                    let cj = faces[i].children[j] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[cj]) {
                        child.f[next(j as i64) as usize] = faces[i].children[3];
                    }

                    // Update children `f` pointers for neighbor children.
                    let f2 = faces[i].f[j];
                    if let Some(child) = Arc::get_mut(&mut new_faces[cj]) {
                        child.f[j] = if f2 != -1 {
                            let face2 = Arc::clone(&faces[f2 as usize]);
                            face2.children[face2.vnum(faces[i].v[j]) as usize]
                        } else {
                            -1
                        };
                    }

                    let f2 = faces[i].f[prev(j as i64) as usize];
                    let cj = faces[i].children[j] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[cj]) {
                        child.f[prev(j as i64) as usize] = if f2 != -1 {
                            let face2 = Arc::clone(&faces[f2 as usize]);
                            face2.children[face2.vnum(faces[i].v[j]) as usize]
                        } else {
                            -1
                        };
                    }
                }
            }

            // Update face vertex pointers
            for face in faces.iter() {
                for j in 0..3_usize {
                    // Update child vertex pointer to new even vertex.
                    let cj = face.children[j] as usize;
                    let mut new_vi = -1_i64;
                    if let Some(child) = Arc::get_mut(&mut new_faces[cj]) {
                        child.v[j] = verts[face.v[j] as usize].child;

                        // Update child vertex pointer to new odd vertex
                        let edge = SDEdge::new(face.v[j], face.v[next(j as i64) as usize]);
                        if let Some(&vi) = edge_verts.get(&edge) {
                            new_vi = vi;
                            child.v[next(j as i64) as usize] = new_vi;
                        }
                    }

                    let cnextj = face.children[next(j as i64) as usize] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[cnextj]) {
                        child.v[j] = new_vi;
                    }

                    let c3 = face.children[3] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[c3]) {
                        child.v[j] = new_vi;
                    }
                }
            }

            // Prepare for next level of subdivision.
            faces = new_faces.split_off(0);
            verts = new_vertices.split_off(0);
        }

        // Push vertices to limit surface.
        let mut p_limit: Vec<Point3f> = Vec::with_capacity(verts.len());
        for i in 0..verts.len() {
            if verts[i].boundary {
                p_limit.push(weight_boundary(
                    Arc::clone(&verts[i]),
                    1.0 / 5.0,
                    i as i64,
                    &verts,
                    &faces,
                ));
            } else {
                p_limit.push(weight_one_ring(
                    Arc::clone(&verts[i]),
                    loop_gamma(Arc::clone(&verts[i]).valence(i as i64, &faces)),
                    i as i64,
                    &verts,
                    &faces,
                ));
            }
        }
        for i in 0..verts.len() {
            Arc::get_mut(&mut verts[i]).unwrap().p = p_limit[i];
        }

        // Compute vertex tangents on limit surface.
        let mut ns: Vec<Normal3f> = Vec::with_capacity(verts.len());
        for i in 0..verts.len() {
            let vertex = Arc::clone(&verts[i]);
            let mut s = Vector3f::default();
            let mut t = Vector3f::default();
            let valence = vertex.valence(i as i64, &faces);
            let p_ring = vertex.one_ring(i as i64, &verts, &faces);
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
                        -1.0 * p_ring[0] + 2.0 * p_ring[1] + 2.0 * p_ring[2]
                            - 1.0 * p_ring[3]
                            - 2.0 * Vector3f::from(vertex.p),
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
        let n_tris = faces.len();
        let mut vertex_indices: Vec<usize> = Vec::with_capacity(3 * n_tris);
        for face in faces.iter() {
            for j in 0..3 {
                vertex_indices.push(face.v[j] as usize);
            }
        }

        TriangleMesh::create(
            Arc::clone(&object_to_world),
            Arc::clone(&world_to_object),
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
            Arc::clone(&o2w),
            Arc::clone(&w2o),
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
/// * `vi`    - The vertex index.
/// * `verts` - The list of vertices.
/// * `faces` - The list of faces.
fn weight_one_ring(
    vert: Arc<SDVertex>,
    beta: Float,
    vi: i64,
    verts: &[Arc<SDVertex>],
    faces: &[Arc<SDFace>],
) -> Point3f {
    let valence = vert.valence(vi, faces);
    let p_ring = vert.one_ring(vi, verts, faces);
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
/// * `vi`    - The vertex index.
/// * `verts` - The list of vertices.
/// * `faces` - The list of faces.
fn weight_boundary(
    vert: Arc<SDVertex>,
    beta: Float,
    vi: i64,
    verts: &[Arc<SDVertex>],
    faces: &[Arc<SDFace>],
) -> Point3f {
    // put _vert_ one-ring in _p_ring_
    let valence = vert.valence(vi, faces);
    let p_ring = vert.one_ring(vi, verts, faces);
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
