//! PLY mesh.

#![allow(dead_code)]

use super::TriangleMesh;
use core::geometry::*;
use core::paramset::*;
use core::texture::FloatTextureMap;
use ply_rs::parser::Parser;
use ply_rs::ply::*;
use std::fs::File;
use std::io::BufReader;
use std::sync::Arc;
use textures::ConstantTexture;

/// Implements PLY mesh loading.
pub struct PLYMesh;

impl PLYMesh {
    /// Create mesh from a PLY file using given parameter set, object to world
    /// transform, world to object transform, whether or not surface normal
    /// orientation is reversed and floating point texture maps.
    ///
    /// NOTE: Because we return a set of shapes as `Vec<Arc<Shape>>` we cannot
    /// implement this as `From` trait :(
    ///
    /// * `p` - A tuple containing the parameter set, object to world transform,
    ///         world to object transform, whether or not surface normal
    ///         orientation is reversed and floating point texture maps.
    pub fn from_props(
        p: (
            &ParamSet,
            ArcTransform,
            ArcTransform,
            bool,
            &FloatTextureMap,
        ),
    ) -> Vec<ArcShape> {
        let (params, o2w, w2o, reverse_orientation, float_textures) = p;

        let path = params.find_one_filename("filename", String::from(""));
        assert!(path.len() > 0, "PLY filename not provied");

        // Look up an alpha texture, if applicable.
        let alpha_tex_name = params.find_one_texture("alpha", String::from(""));
        let alpha_tex = if alpha_tex_name.len() > 0 {
            if let Some(tex) = float_textures.get(&alpha_tex_name) {
                Arc::clone(&tex)
            } else {
                warn!(
                    "Couldn't find float texture '{}' for 'alpha' parameter",
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
                Arc::clone(tex)
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

        // Parse the PLY file.
        let file = File::open(&path).expect(format!("Unable to open PLY file '{}'", path).as_ref());
        let mut reader = BufReader::new(file);

        let parser = Parser::<DefaultElement>::new();
        let ply = match parser.read_ply(&mut reader) {
            Ok(p) => p,
            Err(e) => panic!("Unable to parse PLY file '{}'. {}.", path, e),
        };

        let mut points: Vec<Point3f> = vec![];
        let mut normals: Vec<Normal3f> = vec![];
        let mut uvs: Vec<Point2f> = vec![];
        let mut has_normals = true;
        let mut has_uvs = true;
        let mut vertex_indices: Vec<usize> = vec![];
        let mut face_count = 0;

        for (name, list) in ply.payload.iter() {
            match name.as_ref() {
                "vertex" => {
                    for elem in list.iter() {
                        // Parse the vertex.
                        let vertex = Self::parse_vertex(elem);

                        // Record the vertex point.
                        points.push(vertex.point);

                        // Continue to record normals unless we find one missing.
                        has_normals = has_normals && vertex.has_normal;
                        if has_normals {
                            normals.push(vertex.normal);
                        }

                        // Continue to record normals unless we find one missing.
                        has_uvs = has_uvs && vertex.has_uv;
                        if has_uvs {
                            uvs.push(vertex.uv);
                        }
                    }
                }
                "face" => {
                    for elem in list.iter() {
                        // Parse face.
                        Self::parse_face(elem, &mut vertex_indices);
                        face_count += 1;
                    }
                }
                s => warn!("Ignoring unexpected element '{}' in '{}'", s, path),
            }
        }

        // Inspect the structure of the PLY file.
        if points.len() == 0 || face_count == 0 {
            error!(
                "PLY file '{}' is invalid! No face/vertex elements found!",
                path
            );
            return vec![];
        }

        // If there are missing vertex normals, ignore them entirely.
        if !has_normals {
            normals = vec![];
        }

        // If there are missing UVs, ignore them entirely.
        if !has_uvs {
            uvs = vec![];
        }

        // Create the triangle mesh.
        TriangleMesh::create(
            Arc::clone(&o2w),
            Arc::clone(&w2o),
            reverse_orientation,
            vertex_indices.len() / 3,
            vertex_indices,
            points,
            normals,
            vec![],
            uvs,
            Some(alpha_tex),
            Some(shadow_alpha_tex),
            vec![],
        )
    }

    /// Parse vertex data. Because of the way the parser works it gives one
    /// property at a time. Vertex normals and UV coordinates are treated as
    /// optional. We are assuming if any one property of a vertex, normal or UV
    /// are present then all properties of corresponding data are present. So,
    /// bad data where a property is missing will result in bad return values.
    ///
    /// * `elem` - A map of property names and values.
    fn parse_vertex(elem: &KeyMap<Property>) -> Vertex {
        let mut p = Point3f::default();
        let mut n = Normal3f::default();
        let mut uv = Point2f::default();
        let mut nc = 0;
        let mut uvc = 0;

        for (name, value) in elem.iter() {
            if let Property::Float(v) = value {
                match name.as_ref() {
                    "x" => {
                        p.x = *v;
                    }
                    "y" => {
                        p.y = *v;
                    }
                    "z" => {
                        p.z = *v;
                    }
                    "nx" => {
                        n.x = *v;
                        nc += 1;
                    }
                    "ny" => {
                        n.y = *v;
                        nc += 1;
                    }
                    "nz" => {
                        n.z = *v;
                        nc += 1;
                    }
                    "u" | "s" | "texture_u" | "texture_s" => {
                        uv.x = *v;
                        uvc += 1;
                    }
                    "v" | "t" | "texture_v" | "texture_t" => {
                        uv.y = *v;
                        uvc += 1;
                    }
                    s => debug!("Ignoring unexpected vertex element '{}'", s),
                }
            } else {
                debug!("Ignoring unexpected vertex property type");
            }
        }

        Vertex::new(p, n, uv, nc == 3, uvc == 2)
    }

    /// Parse face data. Only vertex indices are supporeted. Unfortunately,
    /// these can be present as i32 or u32 and have to be converted to usize.
    ///
    /// * `elem` - A map of property names and values.
    fn parse_face(elem: &KeyMap<Property>, vertex_indices: &mut Vec<usize>) {
        for (name, value) in elem.iter() {
            match name.as_ref() {
                "vertex_indices" => {
                    if let Property::ListInt(vi) = value {
                        if vi.len() != 3 && vi.len() != 4 {
                            panic!("Only triangles and quads are supported!");
                        }
                        if vi.len() >= 3 {
                            vertex_indices.push(vi[0] as usize);
                            vertex_indices.push(vi[1] as usize);
                            vertex_indices.push(vi[2] as usize);
                        }
                        if vi.len() == 4 {
                            vertex_indices.push(vi[3] as usize);
                            vertex_indices.push(vi[0] as usize);
                            vertex_indices.push(vi[2] as usize);
                        }
                    } else if let Property::ListUInt(vi) = value {
                        if vi.len() != 3 && vi.len() != 4 {
                            panic!("Only triangles and quads are supported!");
                        }
                        if vi.len() >= 3 {
                            vertex_indices.push(vi[0] as usize);
                            vertex_indices.push(vi[1] as usize);
                            vertex_indices.push(vi[2] as usize);
                        }
                        if vi.len() == 4 {
                            vertex_indices.push(vi[3] as usize);
                            vertex_indices.push(vi[0] as usize);
                            vertex_indices.push(vi[2] as usize);
                        }
                    } else {
                        debug!("Ignoring unexpected face property type");
                    }
                }
                s => debug!("Ignoring unexpected face element '{}'", s),
            }
        }
    }
}

struct Vertex {
    point: Point3f,
    normal: Normal3f,
    uv: Point2f,
    has_normal: bool,
    has_uv: bool,
}

impl Vertex {
    fn new(point: Point3f, normal: Normal3f, uv: Point2f, has_normal: bool, has_uv: bool) -> Self {
        Self {
            point,
            normal,
            uv,
            has_normal,
            has_uv,
        }
    }
}
