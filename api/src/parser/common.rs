//! Common

use crate::{Api, ParamSet};
use core::geometry::{Normal3f, Point2f, Point3f, Vector2f, Vector3f};
use core::pbrt::{Float, Int};
use std::fmt;
use std::fs;

/// Represents the `pbrt` rule of the PBRT file.
pub(crate) struct Pbrt {
    /// Absolute path to scene file.
    pub(crate) scene_path: String,

    /// Parsed `stmt`s.
    pub(crate) stmts: Vec<Stmt>,
}

impl Pbrt {
    /// Process the rule.
    ///
    /// * `api` - The PBRT API interface.
    pub(crate) fn process(&self, api: &mut Api) {
        api.set_current_working_dir(&self.scene_path);
        self.stmts.iter().for_each(|stmt| stmt.process(api));
    }
}

/// Represents the `stmt` rule of the PBRT file.
pub(crate) enum Stmt {
    Include(String),    // include_stmt(abs_path)
    Block(BlockStmt),   // block_stmt
    Option(OptionStmt), // option_stmt
    Scene(SceneStmt),   // scene_stmt
    CTM(CTMStmt),       // ctm_stmt
}

impl Stmt {
    /// Process the rule.
    ///
    /// * `api` - The PBRT API interface.
    pub(crate) fn process(&self, api: &mut Api) {
        match self {
            Self::Include(ref abs_path) => match super::parse(abs_path, api) {
                Ok(()) => (),
                Err(e) => error!("Error parsing include file '{}'. {}", abs_path, e),
            },
            Self::Block(stmt) => stmt.process(api),
            Self::Option(stmt) => stmt.process(api),
            Self::Scene(stmt) => stmt.process(api),
            Self::CTM(stmt) => stmt.process(api),
        }
    }
}

/// Represents the `block_stmt` rule of the PBRT file.
pub(crate) enum BlockStmt {
    WorldBegin,          // world_begin_stmt
    WorldEnd,            // world_end_stmt
    AttributeBegin,      // attribute_begin_stmt
    AttributeEnd,        // attribute_end_stmt
    ObjectBegin(String), // object_begin_stmt(name)
    ObjectEnd,           // object_end_stmt
    TransformBegin,      // transform_begin_stmt
    TransformEnd,        // transform_end_stmt
}

impl BlockStmt {
    /// Process the rule.
    ///
    /// * `api` - The PBRT API interface.
    pub(crate) fn process(&self, api: &mut Api) {
        match self {
            Self::WorldBegin => api.pbrt_world_begin(),
            Self::WorldEnd => api.pbrt_world_end(),
            Self::AttributeBegin => api.pbrt_attribute_begin(),
            Self::AttributeEnd => api.pbrt_attribute_end(),
            Self::ObjectBegin(name) => api.pbrt_object_begin(name.clone()),
            Self::ObjectEnd => api.pbrt_object_end(),
            Self::TransformBegin => api.pbrt_transform_begin(),
            Self::TransformEnd => api.pbrt_transform_end(),
        }
    }
}

/// Represents the `option_stmt` rule of the PBRT file.
pub(crate) enum OptionStmt {
    Accelerator(String, Vec<Param>),     // accelerator_stmt(name, params)
    Camera(String, Vec<Param>),          // camera_stmt(name, params)
    Film(String, Vec<Param>),            // film_stmt(name, params)
    Integrator(String, Vec<Param>),      // integrator_stmt(name, params)
    MakeNamedMedium(String, Vec<Param>), // make_named_medium_stmt(name, params)
    Sampler(String, Vec<Param>),         // sampler_stmt(name, params)
    PixelFilter(String, Vec<Param>),     // pixel_filter_stmt(name, params)
}

impl OptionStmt {
    /// Process the rule.
    ///
    /// * `api` - The PBRT API interface.
    pub(crate) fn process(&self, api: &mut Api) {
        match self {
            Self::Accelerator(name, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_accelerator(name.to_owned(), &ps);
            }
            Self::Camera(name, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_camera(name.to_owned(), &ps);
            }
            Self::Film(name, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_film(name.to_owned(), &ps);
            }
            Self::Integrator(name, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_integrator(name.to_owned(), &ps);
            }
            Self::MakeNamedMedium(name, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_make_named_medium(name.to_owned(), &ps);
            }
            Self::Sampler(name, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_sampler(name.to_owned(), &ps);
            }
            Self::PixelFilter(name, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_pixel_filter(name.to_owned(), &ps);
            }
        }
    }
}

/// Represents the `scene_stmt` rule of the PBRT file.
pub(crate) enum SceneStmt {
    AreaLightSource(String, Vec<Param>), // area_light_source_stmt(name, params)
    LightSource(String, Vec<Param>),     // light_source_stmt(name, params)
    MakeNamedMaterial(String, Vec<Param>), // make_named_material_stmt(name, params)
    Material(String, Vec<Param>),        // material_stmt(name, params)
    Shape(String, Vec<Param>),           // shape_stmt(name, params)
    NamedMaterial(String),               // named_material_stmt(name)
    ObjectInstance(String),              // object_instance_stmt(name)
    Texture(String, String, String, Vec<Param>), // texture_stmt(name, type, class, params)
    ReverseOrientation,                  // reverse_orientation_stmt
    MediumInterface(String, String),     // medium_interface_stmt(inside, outside)
    ActiveTransform(ActiveTransformTime), // active_transform_stmt(time)
}

impl SceneStmt {
    /// Process the rule.
    ///
    /// * `api` - The PBRT API interface.
    pub(crate) fn process(&self, api: &mut Api) {
        match self {
            Self::AreaLightSource(name, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_area_light_source(name.to_owned(), &ps);
            }
            Self::LightSource(name, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_light_source(name.to_owned(), &ps);
            }
            Self::MakeNamedMaterial(name, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_make_named_material(name.to_owned(), &ps);
            }
            Self::Material(name, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_material(name.to_owned(), &ps);
            }
            Self::Shape(name, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_shape(name.to_owned(), &ps);
            }
            Self::NamedMaterial(name) => api.pbrt_named_material(name.to_owned()),
            Self::ObjectInstance(name) => api.pbrt_object_instance(name.to_owned()),
            Self::Texture(name, tex_type, tex_class, params) => {
                let ps = Param::params_to_paramset(params, &api.cwd);
                api.pbrt_texture(
                    name.to_owned(),
                    tex_type.to_owned(),
                    tex_class.to_owned(),
                    &ps,
                );
            }
            Self::ReverseOrientation => api.pbrt_reverse_orientation(),
            Self::MediumInterface(inside, outside) => {
                api.pbrt_medium_interface(inside.to_owned(), outside.to_owned())
            }
            Self::ActiveTransform(time) => match time {
                ActiveTransformTime::Start => api.pbrt_active_transform_start_time(),
                ActiveTransformTime::End => api.pbrt_active_transform_end_time(),
                ActiveTransformTime::All => api.pbrt_active_transform_all(),
            },
        }
    }
}

/// Represents the `ctm_stmt` rule of the PBRT file.
pub(crate) enum CTMStmt {
    Identity,                                   // identity_stmt
    Translate([Float; 3]),                      // translate_stmt(x, y, z)
    Scale([Float; 3]),                          // scale_stmt(x, y, z)
    Rotate([Float; 4]),                         // rotate_stmt(angle, x, y, z)
    LookAt([Float; 3], [Float; 3], [Float; 3]), // look_at_stmt(eye, look_at, up)
    CoordinateSystem(String),                   // coordinate_system_stmt(name)
    CoordinateSystemTransform(String),          // coord_sys_transform_stmt(name)
    Transform([Float; 16]),                     // transform_stmt(matrix)
    ConcatTransform([Float; 16]),               // concat_transform_stmt(matrix)
    TransformTimes([Float; 2]),                 // transform_times(start, end)
}

impl CTMStmt {
    /// Process the rule.
    ///
    /// * `api` - The PBRT API interface.
    pub(crate) fn process(&self, api: &mut Api) {
        match self {
            Self::Identity => api.pbrt_identity(),
            Self::Translate(v) => api.pbrt_translate(v[0], v[1], v[2]),
            Self::Scale(v) => api.pbrt_scale(v[0], v[1], v[2]),
            Self::Rotate(v) => api.pbrt_rotate(v[0], v[1], v[2], v[3]),
            Self::LookAt(eye, look_at, up) => api.pbrt_look_at(
                eye[0], eye[1], eye[2], look_at[0], look_at[1], look_at[2], up[0], up[1], up[2],
            ),
            Self::CoordinateSystem(name) => api.pbrt_coordinate_system(name.clone()),
            Self::CoordinateSystemTransform(name) => api.pbrt_coord_sys_transform(name.clone()),
            Self::Transform(tr) => api.pbrt_transform(&tr),
            Self::ConcatTransform(tr) => api.pbrt_concat_transform(&tr),
            Self::TransformTimes(v) => api.pbrt_transform_times(v[0], v[1]),
        }
    }
}

/// Represents the `active_transform_time` rule of the PBRT file.
pub(crate) enum ActiveTransformTime {
    Start, // StartTime
    End,   // EndTime
    All,   // All
}

/// Represents the `colour_type` rule of the PBRT file.
pub(crate) enum ColourType {
    RGB, // colour, color, rgb
    XYZ, // xyz
}

/// Represents the `param` rule of the PBRT file.
pub(crate) enum Param {
    Point3d(String, Vec<Point3f>),          // point_3d_param(name, p)
    Vector3d(String, Vec<Vector3f>),        // vector_3d_param(name, v)
    Normal3d(String, Vec<Normal3f>),        // normal_3d_param(name, n)
    Point2d(String, Vec<Point2f>),          // point_2d_param(name, p)
    Vector2d(String, Vec<Vector2f>),        // vector_2d_param(name, v)
    String(String, Vec<String>),            // str_param(name, [s])
    Bool(String, Vec<bool>),                // bool_param(name, [b])
    Float(String, Vec<Float>),              // float_param(name, [f])
    Int(String, Vec<Int>),                  // int_param(name, [i])
    Colour(String, ColourType, Vec<Float>), // colour_param(name, type, [c])
    Blackbody(String, Vec<Float>),          // blackbody_param(name, [f])
    SpectrumFloat(String, Vec<Float>),      // spectrum_param(name, [f])
    SpectrumFile(String, String),           // spectrum_param(name, path)
    Texture(String, Vec<String>),           // texture_param(name, [s])
}

impl Param {
    /// Builds a `ParamSet` from a `Vec<Param<'i>>`.
    ///
    /// * `params` - List of params.
    /// * `cwd`    - Current working directory for handling relative paths.
    pub(crate) fn params_to_paramset(params: &[Param], cwd: &str) -> ParamSet {
        let mut ps = ParamSet::new();
        for param in params.iter() {
            match param {
                Self::Point3d(name, v) => ps.add_point3f(&name, &v),
                Self::Vector3d(name, v) => ps.add_vector3f(&name, &v),
                Self::Normal3d(name, v) => ps.add_normal3f(&name, &v),
                Self::Point2d(name, v) => ps.add_point2f(&name, &v),
                Self::Vector2d(name, v) => ps.add_vector2f(&name, &v),
                Self::String(name, v) => ps.add_string(&name, &v),
                Self::Bool(name, v) => ps.add_bool(&name, &v),
                Self::Float(name, v) => ps.add_float(&name, &v),
                Self::Int(name, v) => ps.add_int(&name, &v),
                Self::Colour(name, t, v) => match t {
                    ColourType::RGB => ps.add_rgb_spectrum(&name, &v),
                    ColourType::XYZ => ps.add_xyz_spectrum(&name, &v),
                },
                Self::Blackbody(name, v) => ps.add_blackbody_spectrum(&name, &v),
                Self::SpectrumFloat(name, v) => ps.add_sampled_spectrum(&name, &v),
                Self::SpectrumFile(name, path) => {
                    ps.add_sampled_spectrum_files(&name, &[path.clone()], cwd)
                }
                Self::Texture(name, v) => ps.add_texture(&name, &v),
            }
        }
        ps
    }
}

/// Error parsing of `active_transform_type` rule.
#[derive(Debug, Clone)]
pub(crate) struct ActiveTransformTypeError<'i>(pub(crate) &'i str);

impl<'i> fmt::Display for ActiveTransformTypeError<'i> {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - The formatter.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "invalid active transform type '{}', Allowed values are 'StartTime', 'EndTime', 'All'.",
            self.0
        )
    }
}

/// Error parsing of `colour_type` rule.
#[derive(Debug, Clone)]
pub(crate) struct ColourTypeError<'i>(pub(crate) &'i str);

impl<'i> fmt::Display for ColourTypeError<'i> {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - The formatter.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "invalid colour type '{}', Allowed values are 'colour', 'color', 'rgb', 'xyz'.",
            self.0
        )
    }
}

/// Read the entire file and return its contents as a String.
///
/// * `abs_path` - Absolute path to file.
pub(crate) fn file_to_string(abs_path: &str) -> std::result::Result<String, String> {
    match fs::read_to_string(&abs_path) {
        Ok(s) => Ok(s),
        Err(e) => Err(format!("Error reading file '{}': {}", abs_path, e)),
    }
}
