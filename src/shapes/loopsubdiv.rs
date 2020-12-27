//! Loop Subdivision Surfaces.

#![allow(dead_code)]
use super::TriangleMesh;
use crate::core::geometry::*;
use crate::core::pbrt::*;
use std::cmp::{Ordering, PartialOrd};
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::sync::Arc;

/// Subdivision surface vertex.
#[derive(Clone, Default)]
struct SDVertex {
    /// Index in vertices.
    pub id: usize,

    /// Position.
    pub p: Point3f,

    /// An arbitrary face adjacent to this vertex.
    pub start_face: Option<usize>,

    /// First vertex in next level of subdivision.
    pub child: Option<usize>,

    /// Indicates regular or extraordinary vertex. Interior vertices with valence
    /// other than six, or boundary vertices with valence other than four, are
    /// called extraordinary vertices; otherwise, they are called regular.
    pub regular: bool,

    /// Indicates vertex is on the boundary of the mesh.
    pub boundary: bool,
}

impl SDVertex {
    /// Create a new subdivision vertex.
    fn new(id: usize, p: &Point3f) -> Self {
        Self {
            id,
            p: *p,
            start_face: None,
            child: None,
            regular: false,
            boundary: false,
        }
    }

    /// Returns number of vertices directly adjacent to this vertex.
    pub fn valence(&self, faces: &[Arc<SDFace>]) -> usize {
        if !self.boundary {
            // Compute valence of interior vertex.
            let mut face = self.start_face;
            let mut nf = 1_usize;
            while let Some(f) = face {
                face = faces[f].next_face(self.id);
                if face.is_none() || face == self.start_face {
                    break;
                }
                nf += 1;
            }
            nf
        } else {
            // Compute valence of boundary vertex.
            let mut face = self.start_face;
            let mut nf = 1_usize;
            while let Some(f) = face {
                face = faces[f].next_face(self.id);
                if face.is_none() || face == self.start_face {
                    break;
                }
                nf += 1;
            }
            face = self.start_face;
            while let Some(f) = face {
                face = faces[f].prev_face(self.id);
                if face.is_none() || face == self.start_face {
                    break;
                }
                nf += 1;
            }
            nf + 1
        }
    }

    /// Returns position of vertices around this vertex.
    ///
    /// * `verts` - The vertices.
    /// * `faces` - The faces.
    pub fn one_ring(&self, verts: &[Arc<SDVertex>], faces: &[Arc<SDFace>]) -> Vec<Point3f> {
        let mut p: Vec<Point3f> = vec![];

        if !self.boundary {
            // Get one-ring vertices from interior vertex.
            let mut face = self.start_face;
            while let Some(f) = face {
                let next_vert = faces[f].next_vert(self.id);
                p.push(verts[next_vert].p);

                face = faces[f].next_face(self.id);
                if face == self.start_face {
                    break;
                }
            }
        } else {
            // Get one-ring vertices from boundary vertex.
            let mut face = self.start_face;
            while let Some(f) = face {
                if let Some(face2) = faces[f].next_face(self.id) {
                    face = Some(face2);
                } else {
                    break;
                }
            }

            if let Some(mut f) = face {
                p.push(verts[faces[f].next_vert(self.id)].p);
                loop {
                    let prev_vert = faces[f].prev_vert(self.id);
                    p.push(verts[prev_vert].p);

                    if let Some(f2) = faces[f].prev_face(self.id) {
                        f = f2;
                    } else {
                        break;
                    }
                }
            }
        }

        p
    }
}

impl PartialEq for SDVertex {
    /// Determines equality based on the position.
    fn eq(&self, other: &Self) -> bool {
        self.p == other.p
    }
}

/// Subdivision surface edge.
#[derive(Clone, Default, Eq)]
struct SDEdge {
    /// Vertex indices for the endpoints of the edge.
    pub v: [usize; 2],

    /// The face indices of adjacent faces.
    pub f: [Option<usize>; 2],

    /// The index of first face found that is adjacent to this edge.
    pub f0_edge_num: Option<usize>,
}

impl SDEdge {
    /// Create a new subdvision surface edge.
    ///
    /// * `v0` - First endpoint.
    /// * `v1` - Second endpoint.
    fn new(v0: usize, v1: usize) -> Self {
        Self {
            v: [v0, v1],
            f: [None, None],
            f0_edge_num: None,
        }
    }
}

impl PartialEq for SDEdge {
    /// Determine equality.
    fn eq(&self, other: &Self) -> bool {
        self.v[0] == other.v[0]
    }
}

impl PartialOrd for SDEdge {
    /// Define an ordering for edges.
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
    // Feeds this value into the given `Hasher`.
    //
    // * `state` - Hasher state.
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.v[0].hash(state);
        self.v[1].hash(state);
    }
}

/// Returns the next index modulo 3.
macro_rules! next {
    ($i: expr) => {
        ((($i) + 1) % 3)
    };
}

/// Returns the previous index modulo 3.
macro_rules! prev {
    ($i: expr) => {
        ((($i) + 2) % 3)
    };
}

/// Subdivision surface face.
#[derive(Clone, Default)]
struct SDFace {
    /// Index of vertices of the face.
    pub v: [usize; 3],

    /// Adjacent face indices.
    pub f: [Option<usize>; 3],

    /// Indices of faces at the next level of subdivision.
    pub children: [usize; 4],
}

impl SDFace {
    /// Create a new subdivision face from vertex indices.
    ///
    /// * `v`: A slice containing 3 vertex indices.
    fn new(v: &[usize]) -> Self {
        debug_assert!(v.len() == 3, "invalid slice of length != 3");
        Self {
            v: [v[0], v[1], v[2]],
            f: [None, None, None],
            children: [0, 0, 0, 0],
        }
    }

    /// Find the index in v[] of a given vertex.
    ///
    /// * `vert` - The vertex index.
    pub fn vnum(&self, vert: usize) -> usize {
        self.v
            .iter()
            .position(|&v| v == vert)
            .expect("Basic logic error in SDFace::vnum()")
    }

    /// Returns index of next face adjacent to given vertex.
    ///
    /// * `vert` - The vertex index.
    pub fn next_face(&self, vert: usize) -> Option<usize> {
        self.f[self.vnum(vert)]
    }

    /// Returns index of previous face adjacent to given vertex.
    ///
    /// * `vert` - The vertex index.
    pub fn prev_face(&self, vert: usize) -> Option<usize> {
        self.f[prev!(self.vnum(vert))]
    }

    /// Returns index of next vertex adjacent to given vertex.
    ///
    /// * `vert` - The vertex index.
    pub fn next_vert(&self, vert: usize) -> usize {
        self.v[next!(self.vnum(vert))]
    }

    /// Returns index of previous vertex adjacent to given vertex.
    ///
    /// * `vert` - The vertex index.
    pub fn prev_vert(&self, vert: usize) -> usize {
        self.v[prev!(self.vnum(vert))]
    }

    /// Returns the vertex index of vertex opposite an edge.
    ///
    /// * `v0` - Vertex index of first endpoint of edge.
    /// * `v1` - Vertex index of second endpoint of edge.
    pub fn other_vert(&self, v0: usize, v1: usize) -> usize {
        *self
            .v
            .iter()
            .find(|&&v| v != v0 && v != v1)
            .expect("Basic logic error in SDVertex::other_vert()")
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
        // Create vertices from points.
        let mut vertices: Vec<Arc<SDVertex>> = p
            .iter()
            .enumerate()
            .map(|(id, point)| Arc::new(SDVertex::new(id, &point)))
            .collect();

        // Create faces and set vertices.
        let n_faces = vertex_indices.len() / 3;
        let mut faces: Vec<Arc<SDFace>> = (0..n_faces)
            .map(|i| {
                let fv: Vec<usize> = (0..3)
                    .map(|j| i * 3 * j)
                    .map(|vert_idx| {
                        Arc::get_mut(&mut vertices[vert_idx]).unwrap().start_face = Some(i);
                        vert_idx
                    })
                    .collect();
                Arc::new(SDFace::new(&fv))
            })
            .collect();

        // Set neighbour pointers in faces.
        let mut edges: HashSet<SDEdge> = HashSet::new();
        for i in 0..n_faces {
            for edge_num in 0..3 {
                let (v0, v1) = (edge_num, next!(edge_num));

                let mut e = SDEdge::new(faces[i].v[v0], faces[i].v[v1]);

                if let Some(e1) = edges.take(&e) {
                    // Handle previously seen edge. The `take()` will remove it
                    // from the set as well which we want because we'll have found
                    // the second face on this edge.
                    if let (Some(f0), Some(f0_edge_num)) = (e1.f[0], e1.f0_edge_num) {
                        if let Some(face_f0) = Arc::get_mut(&mut faces[f0]) {
                            face_f0.f[f0_edge_num] = Some(i);
                        } else {
                            panic!("loop_subdivide(): faces[f0={}] not found", f0);
                        }

                        if let Some(face_i) = Arc::get_mut(&mut faces[i]) {
                            face_i.f[edge_num] = e1.f[0];
                        } else {
                            panic!("loop_subdivide(): faces[i={}] not found.", i);
                        }
                    } else {
                        panic!("loop_subdivide(): missing e1.f[0] or e1.f0_edge_num.");
                    }
                } else {
                    // Handle new edge.
                    e.f[0] = Some(i);
                    e.f0_edge_num = Some(edge_num);
                    edges.insert(e);
                }
            }
        }

        // Finish vertex initialization.
        for (i, v) in vertices.iter_mut().enumerate() {
            let mut face: Option<usize> = v.start_face;
            while let Some(f) = face {
                face = faces[f].next_face(i);
                if face == v.start_face {
                    break;
                }
            }

            if let Some(vert) = Arc::get_mut(v) {
                vert.boundary = face.is_none();
                vert.regular = (!vert.boundary && vert.valence(&faces) == 6)
                    || (vert.boundary && vert.valence(&faces) == 4);
            } else {
                panic!(
                    "loop_subdivide(): unable to finish initialization of vertices[i={}].",
                    i
                );
            }
        }

        // Refine loop subdivision into triangles.
        for level in 0..n_levels {
            // Update next level of subdivision.
            let mut new_faces: Vec<Arc<SDFace>> = vec![];
            let mut new_vertices: Vec<Arc<SDVertex>> = vec![];

            // Allocate next level of children in mesh tree.
            for (i, vertex) in vertices.iter_mut().enumerate() {
                if let Some(vert) = Arc::get_mut(vertex) {
                    vert.child = Some(i);
                    new_vertices.push(Arc::new(SDVertex {
                        id: i,
                        p: Point3f::default(),
                        regular: vert.regular,
                        boundary: vert.boundary,
                        start_face: None,
                        child: None,
                    }));
                } else {
                    panic!(
                        "loop_subdivide(): unable to set vertices[i={}].child at level={}.",
                        i, level
                    );
                }
            }

            for face in faces.iter_mut() {
                for k in 0..4 {
                    if let Some(fc) = Arc::get_mut(face) {
                        fc.children[k] = new_faces.len();
                        new_faces.push(Arc::new(SDFace::default()));
                    } else {
                        panic!(
                            "loop_subdivide(): unable to set face.children[k={}] at level={}.",
                            k, level
                        );
                    }
                }
            }

            // Update vertex positions and create new edge vertices.
            // Update vertex positions for even vertices.
            for vertex in vertices.iter() {
                // Child was already set earlier so this shouldn't fail.
                if let Some(child) = vertex.child {
                    if let Some(child_vertex) = Arc::get_mut(&mut new_vertices[child]) {
                        child_vertex.p = if !vertex.boundary {
                            // Apply one-ring rule for even vertex.
                            if vertex.regular {
                                weight_one_ring(&vertex, 1.0 / 16.0, &vertices, &faces)
                            } else {
                                weight_one_ring(
                                    &vertex,
                                    beta(vertex.valence(&faces)),
                                    &vertices,
                                    &faces,
                                )
                            }
                        } else {
                            weight_boundary(&vertex, 1.0 / 8.0, &vertices, &faces)
                        };
                    } else {
                        panic!(
                        "loop_subdivide(): unable to update new_vertices[child={}].p at level={}.",
                        child, level
                    );
                    }
                }
            }

            // Compute new odd edge vertices.
            let mut edge_verts: HashMap<SDEdge, usize> = HashMap::new();
            for face in faces.iter() {
                for k in 0..3 {
                    // Compute odd vertex on k^th edge.
                    let edge = SDEdge::new(face.v[k], face.v[next!(k)]);
                    if !edge_verts.contains_key(&edge) {
                        let id = new_vertices.len();
                        let boundary = face.f[k].is_none();

                        // Apply edge rules to compute new vertex position.
                        let p = if boundary {
                            0.5 * vertices[edge.v[0]].p + 0.5 * vertices[edge.v[1]].p
                        } else if let Some(fk) = face.f[k] {
                            let face_k = &faces[fk];
                            3.0 / 8.0 * vertices[edge.v[0]].p
                                + 3.0 / 8.0 * vertices[edge.v[1]].p
                                + 1.0 / 8.0 * vertices[face.other_vert(edge.v[0], edge.v[1])].p
                                + 1.0 / 8.0 * vertices[face_k.other_vert(edge.v[0], edge.v[1])].p
                        } else {
                            panic!(
                                "loop_subdivide(): unable to get face.f[k={}] at level={}.",
                                k, level
                            );
                        };

                        // Create and initialize new odd vertex.
                        new_vertices.push(Arc::new(SDVertex {
                            id,
                            p,
                            regular: true,
                            boundary,
                            start_face: Some(face.children[3]),
                            child: None,
                        }));

                        edge_verts.insert(edge, id);
                    }
                }
            }

            // Update new mesh topology.
            // Update even vertex face pointers.
            for vertex in vertices.iter() {
                if let Some(start_face) = vertex.start_face.map(|f| faces[f].clone()) {
                    let vert_num = start_face.vnum(vertex.id);
                    if let Some(child) = vertex
                        .child
                        .map(|child_id| Arc::get_mut(&mut new_vertices[child_id]))
                        .flatten()
                    {
                        child.start_face = Some(start_face.children[vert_num]);
                    }
                }
            }

            // Update face neighbour pointers.
            for face in faces.iter() {
                for j in 0..3 {
                    // Update children f pointers for siblings.
                    if let Some(child) = Arc::get_mut(&mut new_faces[face.children[3]]) {
                        child.f[j] = Some(face.children[next!(j)]);
                    }
                    if let Some(child) = Arc::get_mut(&mut new_faces[face.children[j]]) {
                        child.f[next!(j)] = Some(face.children[3]);
                    }

                    // Update children f pointers for neighbour children.
                    if let Some(child) = Arc::get_mut(&mut new_faces[face.children[j]]) {
                        child.f[j] = face.f[j].map(|f2| {
                            let face2 = faces[f2].clone();
                            face2.children[face2.vnum(face.v[j])]
                        });

                        child.f[prev!(j)] = face.f[prev!(j)].map(|f2| {
                            let face2 = faces[f2].clone();
                            face2.children[face2.vnum(face.v[j])]
                        });
                    }
                }
            }

            // Update face vertex pointers.
            for face in faces.iter() {
                for j in 0..3 {
                    // Update child vertex pointer to new even vertex.
                    if let Some(child) = Arc::get_mut(&mut new_faces[face.children[j]]) {
                        if let Some(v) = vertices[face.v[j]].child {
                            child.v[j] = v;
                        }
                    }

                    // Update child vertex pointer to new odd vertex
                    let edge = SDEdge::new(face.v[j], face.v[next!(j)]);
                    if let Some(vert) = edge_verts.get(&edge) {
                        if let Some(child) = Arc::get_mut(&mut new_faces[face.children[j]]) {
                            child.v[next!(j)] = *vert;
                        }
                        if let Some(child) = Arc::get_mut(&mut new_faces[face.children[next!(j)]]) {
                            child.v[j] = *vert;
                        }
                        if let Some(child) = Arc::get_mut(&mut new_faces[face.children[3]]) {
                            child.v[j] = *vert;
                        }
                    } else {
                        panic!("loop_subdivide(): edge_verts[edge={:?}] not found", edge.v);
                    }
                }
            }

            // Prepare for next level of subdivision
            faces = new_faces.split_off(0);
            vertices = new_vertices.split_off(0);
        }

        // Push vertices to limit surface.
        let p_limit: Vec<Point3f> = vertices
            .iter()
            .map(|v| {
                if v.boundary {
                    weight_boundary(v, 1.0 / 5.0, &vertices, &faces)
                } else {
                    weight_one_ring(v, loop_gamma(v.valence(&faces)), &vertices, &faces)
                }
            })
            .collect();
        for (i, mut v) in vertices.iter_mut().enumerate() {
            Arc::get_mut(&mut v).unwrap().p = p_limit[i];
        }

        // Compute vertex tangents on limit surface.
        let ns: Vec<Normal3f> = vertices
            .iter()
            .map(|v| {
                let valence = v.valence(&faces);
                let p_ring = v.one_ring(&vertices, &faces);
                let (s, t) = if !v.boundary {
                    // Compute tangents of interior face.
                    (0..valence).fold((Vector3f::default(), Vector3f::default()), |(s, t), j| {
                        let theta = TWO_PI * j as Float / valence as Float;
                        (
                            s + (theta).cos() * Vector3::from(p_ring[j]),
                            t + (theta).sin() * Vector3::from(p_ring[j]),
                        )
                    })
                } else {
                    // Compute tangents of boundary face.
                    let s = p_ring[valence - 1] - p_ring[0];
                    let t = match valence {
                        2 => Vector3::from(p_ring[0] + p_ring[1] - 2.0 * v.p),
                        3 => p_ring[1] - v.p,
                        4 => {
                            // regular
                            Vector3::from(
                                -1.0 * p_ring[0]
                                    + 2.0 * p_ring[1]
                                    + 2.0 * p_ring[2]
                                    + -1.0 * p_ring[3]
                                    + -2.0 * v.p,
                            )
                        }
                        _ => {
                            let theta = PI / (valence - 1) as Float;
                            let t1 = (1..valence - 1).fold(
                                Vector3::from(theta.sin() * (p_ring[0] + p_ring[valence - 1])),
                                |t, k| {
                                    let wt =
                                        (2.0 * theta.cos() - 2.0) * ((k as Float) * theta).sin();
                                    t + Vector3::from(wt * p_ring[k])
                                },
                            );
                            -t1
                        }
                    };

                    (s, t)
                };
                Normal3::from(s.cross(&t))
            })
            .collect();

        // Create triangle mesh from subdivision mesh.
        let tot_verts = vertices.len();
        let mut vertex_indices = Vec::<usize>::with_capacity(tot_verts);
        for face in faces.iter() {
            for j in 0..3 {
                vertex_indices.push(face.v[j]);
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
/// * `vert` - The vertex.
/// * `beta` - The weight.
/// * `verts` - The vertices.
/// * `faces` - The faces.
fn weight_one_ring(
    vert: &SDVertex,
    beta: Float,
    verts: &[Arc<SDVertex>],
    faces: &[Arc<SDFace>],
) -> Point3f {
    // Put vert one-ring in p_ring
    let valence = vert.valence(faces);
    let p_ring = vert.one_ring(verts, faces);
    (0..valence).fold((1.0 - valence as Float * beta) * vert.p, |p, i| {
        p + beta * p_ring[i]
    })
}

/// Applies given weights at the boundary vertex and returns a
/// new position.
///
/// * `vert` - The vertex.
/// * `beta` - The weight.
/// * `verts` - The vertices.
/// * `faces` - The faces.
fn weight_boundary(
    vert: &SDVertex,
    beta: Float,
    verts: &[Arc<SDVertex>],
    faces: &[Arc<SDFace>],
) -> Point3f {
    // Put vert one-ring in p_ring
    let valence = vert.valence(faces);
    let p_ring = vert.one_ring(verts, faces);
    (1.0 - 2.0 * beta) * vert.p + beta * p_ring[0] + beta * p_ring[valence - 1]
}
