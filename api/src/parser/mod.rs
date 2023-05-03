//! PBRT File Parser

mod common;

use super::*;
use common::*;
use core::fileutil::*;
use pest::error::ErrorVariant::CustomError;
use pest_consume::{match_nodes, Error, Parser};

type Result<T> = std::result::Result<T, Error<Rule>>;
type Node<'i> = pest_consume::Node<'i, Rule, String>;

/// The `pest_consume::Parser` generated from a grammar.
#[derive(Parser)]
#[grammar = "parser/grammar.pest"]
struct PbrtParser;

#[pest_consume::parser]
impl PbrtParser {
    /// Parse rule `pbrt`.
    fn pbrt(input: Node) -> Result<Pbrt> {
        Ok(match_nodes!(input.into_children();
            [stmt(stmts).., EOI(_)] => Pbrt { stmts: stmts.collect() },
        ))
    }

    /// Parse rule `EOI (end-of-input)`.
    fn EOI(_input: Node) -> Result<()> {
        Ok(())
    }

    /// Parse rule `stmt`.
    fn stmt(input: Node) -> Result<Stmt> {
        Ok(match_nodes!(input.into_children();
            [include_stmt(stmt)] => stmt,
            [block_stmt(stmt)] => Stmt::Block(stmt),
            [option_stmt(stmt)] => Stmt::Option(stmt),
            [scene_stmt(stmt)] => Stmt::Scene(stmt),
            [ctm_stmt(stmt)] => Stmt::CTM(stmt),
            [comment_stmt(_stmt)] => Stmt::Comment
        ))
    }

    /// Parse rule `comment_stmt`.
    fn comment_stmt(input: Node) -> Result<Stmt> {
        Ok(Stmt::Comment)
    }

    /// Parse rule `include_stmt`.
    fn include_stmt(input: Node) -> Result<Stmt> {
        let scene_path = input.user_data().clone();
        let span = input.as_span();

        let (path, params) = match_nodes!(input.into_children();
            [quoted_str_expr(path)] => (path, "".to_string()),
            [quoted_str_expr(path), params] => {
                // We will be injecting unparsed parameters into the included unparsed file contents.
                let node: Node = params;
                let unparsed_params = node.as_str().to_string();
                (path, unparsed_params)
            },
        );
        if is_relative_path(&path) && !scene_path.is_empty() {
            let p = format!("{}/{}", scene_path, path);
            match absolute_path(&p) {
                Ok(abs_path) => Ok(Stmt::Include(abs_path, params)),
                Err(e) => Err(Error::new_from_span(CustomError { message: e }, span)),
            }
        } else {
            Ok(Stmt::Include(path, params))
        }
    }

    /// Parse rule `block_stmt`.
    fn block_stmt(input: Node) -> Result<BlockStmt> {
        Ok(match_nodes!(input.into_children();
            [world_begin_stmt(stmt)] => stmt,
            [world_end_stmt(stmt)] => stmt,
            [attribute_begin_stmt(stmt)] => stmt,
            [attribute_end_stmt(stmt)] => stmt,
            [object_begin_stmt(stmt)] => stmt,
            [object_end_stmt(stmt)] => stmt,
            [transform_begin_stmt(stmt)] => stmt,
            [transform_end_stmt(stmt)] => stmt,
        ))
    }

    /// Parse rule `world_begin_stmt`.
    fn world_begin_stmt(input: Node) -> Result<BlockStmt> {
        Ok(BlockStmt::WorldBegin)
    }

    /// Parse rule `world_end_stmt`.
    fn world_end_stmt(input: Node) -> Result<BlockStmt> {
        Ok(BlockStmt::WorldEnd)
    }

    /// Parse rule `attribute_begin_stmt`.
    fn attribute_begin_stmt(input: Node) -> Result<BlockStmt> {
        Ok(BlockStmt::AttributeBegin)
    }

    /// Parse rule `attribute_end_stmt`.
    fn attribute_end_stmt(input: Node) -> Result<BlockStmt> {
        Ok(BlockStmt::AttributeEnd)
    }

    /// Parse rule `transform_begin_stmt`.
    fn transform_begin_stmt(input: Node) -> Result<BlockStmt> {
        Ok(BlockStmt::TransformBegin)
    }

    /// Parse rule `transform_end_stmt`.
    fn transform_end_stmt(input: Node) -> Result<BlockStmt> {
        Ok(BlockStmt::TransformEnd)
    }

    /// Parse rule `object_begin_stmt`.
    fn object_begin_stmt(input: Node) -> Result<BlockStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => BlockStmt::ObjectBegin(name),
        ))
    }

    /// Parse rule `object_end_stmt`.
    fn object_end_stmt(input: Node) -> Result<BlockStmt> {
        Ok(BlockStmt::ObjectEnd)
    }

    /// Parse rule `option_stmt`.
    fn option_stmt(input: Node) -> Result<OptionStmt> {
        Ok(match_nodes!(input.into_children();
            [accelerator_stmt(stmt)] => stmt,
            [camera_stmt(stmt)] => stmt,
            [film_stmt(stmt)] => stmt,
            [integrator_stmt(stmt)] => stmt,
            [make_named_medium_stmt(stmt)] => stmt,
            [sampler_stmt(stmt)] => stmt,
            [pixel_filter_stmt(stmt)] => stmt,
        ))
    }

    /// Parse rule `accelerator_stmt`.
    fn accelerator_stmt(input: Node) -> Result<OptionStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => OptionStmt::Accelerator(name, vec![]),
            [quoted_str_expr(name), param_list(params)] => OptionStmt::Accelerator(name, params),
        ))
    }

    /// Parse rule `camera_stmt`.
    fn camera_stmt(input: Node) -> Result<OptionStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => OptionStmt::Camera(name, vec![]),
            [quoted_str_expr(name), param_list(params)] => OptionStmt::Camera(name, params),
        ))
    }

    /// Parse rule `film_stmt`.
    fn film_stmt(input: Node) -> Result<OptionStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => OptionStmt::Film(name, vec![]),
            [quoted_str_expr(name), param_list(params)] => OptionStmt::Film(name, params),
        ))
    }

    /// Parse rule `integrator_stmt`.
    fn integrator_stmt(input: Node) -> Result<OptionStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => OptionStmt::Integrator(name, vec![]),
            [quoted_str_expr(name), param_list(params)] => OptionStmt::Integrator(name, params),
        ))
    }

    /// Parse rule `make_named_medium_stmt`.
    fn make_named_medium_stmt(input: Node) -> Result<OptionStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => OptionStmt::MakeNamedMedium(name, vec![]),
            [quoted_str_expr(name), param_list(params)] => OptionStmt::MakeNamedMedium(name, params),
        ))
    }

    /// Parse rule `sampler_stmt`.
    fn sampler_stmt(input: Node) -> Result<OptionStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => OptionStmt::Sampler(name, vec![]),
            [quoted_str_expr(name), param_list(params)] => OptionStmt::Sampler(name, params),
        ))
    }

    /// Parse rule `pixel_filter_stmt`.
    fn pixel_filter_stmt(input: Node) -> Result<OptionStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => OptionStmt::PixelFilter(name, vec![]),
            [quoted_str_expr(name), param_list(params)] => OptionStmt::PixelFilter(name, params),
        ))
    }

    /// Parse rule `scene_stmt`.
    fn scene_stmt(input: Node) -> Result<SceneStmt> {
        Ok(match_nodes!(input.into_children();
            [area_light_source_stmt(stmt)] => stmt,
            [light_source_stmt(stmt)] => stmt,
            [make_named_material_stmt(stmt)] => stmt,
            [material_stmt(stmt)] => stmt,
            [shape_stmt(stmt)] => stmt,
            [texture_stmt(stmt)] => stmt,
            [named_material_stmt(stmt)] => stmt,
            [object_instance_stmt(stmt)] => stmt,
            [reverse_orientation_stmt(stmt)] => stmt,
            [medium_interface_stmt(stmt)] => stmt,
            [active_transform_stmt(stmt)] => stmt,
        ))
    }

    /// Parse rule `area_light_source_stmt`.
    fn area_light_source_stmt(input: Node) -> Result<SceneStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => SceneStmt::AreaLightSource(name, vec![]),
            [quoted_str_expr(name), param_list(params)] => SceneStmt::AreaLightSource(name, params),
        ))
    }

    /// Parse rule `light_source_stmt`.
    fn light_source_stmt(input: Node) -> Result<SceneStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => SceneStmt::LightSource(name, vec![]),
            [quoted_str_expr(name), param_list(params)] => SceneStmt::LightSource(name, params),
        ))
    }

    /// Parse rule `make_named_material_stmt`.
    fn make_named_material_stmt(input: Node) -> Result<SceneStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => SceneStmt::MakeNamedMaterial(name, vec![]),
            [quoted_str_expr(name), param_list(params)] => SceneStmt::MakeNamedMaterial(name, params),
        ))
    }

    /// Parse rule `material_stmt`.
    fn material_stmt(input: Node) -> Result<SceneStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => SceneStmt::Material(name, vec![]),
            [quoted_str_expr(name), param_list(params)] => SceneStmt::Material(name, params),
        ))
    }

    /// Parse rule `shape_stmt`.
    fn shape_stmt(input: Node) -> Result<SceneStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(name)] => SceneStmt::Shape(name, vec![]),
            [quoted_str_expr(name), param_list(params)] => SceneStmt::Shape(name, params),
        ))
    }

    /// Parse rule `texture_stmt`.
    fn texture_stmt(input: Node) -> Result<SceneStmt> {
        Ok(match_nodes!(input.into_children();
            [
                quoted_str_expr(name),
                quoted_str_expr(tex_type),
                quoted_str_expr(tex_class)
            ] => SceneStmt::Texture(name, tex_type, tex_class, vec![]),
            [
                quoted_str_expr(name),
                quoted_str_expr(tex_type),
                quoted_str_expr(tex_class),
                param_list(params)
            ] => SceneStmt::Texture(name, tex_type, tex_class, params),
        ))
    }

    /// Parse rule `named_material_stmt`.
    fn named_material_stmt(input: Node) -> Result<SceneStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_ident_expr(name)] => SceneStmt::NamedMaterial(name),
        ))
    }

    /// Parse rule `object_instance_stmt`.
    fn object_instance_stmt(input: Node) -> Result<SceneStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_ident_expr(name)] => SceneStmt::ObjectInstance(name),
        ))
    }

    /// Parse rule `reverse_orientation_stmt`.
    fn reverse_orientation_stmt(input: Node) -> Result<SceneStmt> {
        Ok(SceneStmt::ReverseOrientation)
    }

    /// Parse rule `medium_interface_stmt`.
    fn medium_interface_stmt(input: Node) -> Result<SceneStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_str(inside)] => SceneStmt::MediumInterface(inside, "".to_owned()),
            [quoted_str(inside), quoted_str(outside)] => SceneStmt::MediumInterface(inside, outside),
        ))
    }

    /// Parse rule `active_transform_stmt`.
    fn active_transform_stmt(input: Node) -> Result<SceneStmt> {
        Ok(match_nodes!(input.into_children();
            [active_transform_time(time)] => SceneStmt::ActiveTransform(time),
        ))
    }

    /// Parse rule `active_transform_time`.
    fn active_transform_time(input: Node) -> Result<ActiveTransformTime> {
        match input.as_str() {
            "StartTime" => Ok(ActiveTransformTime::Start),
            "EndTime" => Ok(ActiveTransformTime::End),
            "All" => Ok(ActiveTransformTime::All),
            s => Err(input.error(ActiveTransformTypeError(s))),
        }
    }

    /// Parse rule `ctm_stmt`.
    fn ctm_stmt(input: Node) -> Result<CTMStmt> {
        Ok(match_nodes!(input.into_children();
            [identity_stmt(stmt)] => stmt,
            [translate_stmt(stmt)] => stmt,
            [scale_stmt(stmt)] => stmt,
            [rotate_stmt(stmt)] => stmt,
            [look_at_stmt(stmt)] => stmt,
            [coordinate_system_stmt(stmt)] => stmt,
            [coord_sys_transform_stmt(stmt)] => stmt,
            [transform_stmt(stmt)] => stmt,
            [concat_transform_stmt(stmt)] => stmt,
            [transform_times_stmt(stmt)] => stmt,
        ))
    }

    /// Parse rule `identity_stmt`.
    fn identity_stmt(input: Node) -> Result<CTMStmt> {
        Ok(CTMStmt::Identity)
    }

    /// Parse rule `translate_stmt`.
    fn translate_stmt(input: Node) -> Result<CTMStmt> {
        Ok(match_nodes!(input.into_children();
            [float_expr(x), float_expr(y), float_expr(z)] =>
                CTMStmt::Translate([x, y, z]),
        ))
    }

    /// Parse rule `scale_stmt`.
    fn scale_stmt(input: Node) -> Result<CTMStmt> {
        Ok(match_nodes!(input.into_children();
            [float_expr(x), float_expr(y), float_expr(z)] =>
                CTMStmt::Scale([x, y, z]),
        ))
    }

    /// Parse rule `rotate_stmt`.
    fn rotate_stmt(input: Node) -> Result<CTMStmt> {
        Ok(match_nodes!(input.into_children();
            [float_expr(angle), float_expr(x), float_expr(y), float_expr(z)] =>
                CTMStmt::Rotate([angle, x, y, z]),
        ))
    }

    /// Parse rule `look_atstmt`.
    fn look_at_stmt(input: Node) -> Result<CTMStmt> {
        Ok(match_nodes!(input.into_children();
            [
                float_expr(v0), float_expr(v1), float_expr(v2), // Eye
                float_expr(v3), float_expr(v4), float_expr(v5), // Look at point
                float_expr(v6), float_expr(v7), float_expr(v8), // Up vector
            ] => CTMStmt::LookAt([v0, v1, v2], [v3, v4, v5], [v6, v7, v8]),
        ))
    }

    /// Parse rule `coordinate_system_stmt`.
    fn coordinate_system_stmt(input: Node) -> Result<CTMStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_ident_expr(name)] => CTMStmt::CoordinateSystem(name),
        ))
    }

    /// Parse rule `coord_sys_transform_stmt`.
    fn coord_sys_transform_stmt(input: Node) -> Result<CTMStmt> {
        Ok(match_nodes!(input.into_children();
            [quoted_ident_expr(name)] => CTMStmt::CoordinateSystemTransform(name),
        ))
    }

    /// Parse rule `transform_stmt`.
    fn transform_stmt(input: Node) -> Result<CTMStmt> {
        let span = input.as_span();
        match_nodes!(input.into_children();
            [float_list_expr(v)] => {
                if v.len() == 16 {
                    let mut values = [0.0 as Float; 16];
                    v.iter().enumerate().for_each(|(i, val)| {
                        values[i] = *val;
                    });
                    Ok(CTMStmt::Transform(values))
                } else {
                    Err(Error::new_from_span(
                        CustomError {
                            message: "transform_stmt does not have neough values".to_owned()
                        },
                        span
                    ))
                }
            }
        )
    }

    /// Parse rule `concat_transform_stmt`.
    fn concat_transform_stmt(input: Node) -> Result<CTMStmt> {
        let span = input.as_span();
        match_nodes!(input.into_children();
            [float_list_expr(v)] => {
                if v.len() == 16 {
                    let mut values = [0.0 as Float; 16];
                    v.iter().enumerate().for_each(|(i, val)| {
                        values[i] = *val;
                    });
                    Ok(CTMStmt::ConcatTransform(values))
                } else {
                    Err(Error::new_from_span(
                        CustomError {
                            message: "concat_transform_stmt does not have neough values".to_owned()
                        },
                        span
                    ))
                }
            }
        )
    }

    /// Parse rule `transform_times_stmt`.
    fn transform_times_stmt(input: Node) -> Result<CTMStmt> {
        Ok(match_nodes!(input.into_children();
            [float(start), float_expr(end)] => CTMStmt::TransformTimes([start, end]),
        ))
    }

    /// Parse rule `param_list`.
    fn param_list(input: Node) -> Result<Vec<Param>> {
        Ok(match_nodes!(input.into_children();
            [param(params)..] => params.collect(),
        ))
    }

    /// Parse rule `param`.
    fn param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [point3d_param(param)] => param,
            [vector3d_param(param)] => param,
            [normal3d_param(param)] => param,
            [point2d_param(param)] => param,
            [vector2d_param(param)] => param,
            [string_param(param)] => param,
            [bool_param(param)] => param,
            [float_param(param)] => param,
            [int_param(param)] => param,
            [colour_param(param)] => param,
            [blackbody_param(param)] => param,
            [spectrum_param(param)] => param,
            [texture_param(param)] => param,
        ))
    }

    /// Parse rule `point3d_param`.
    fn point3d_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [ident(name), float_list_expr(v)] => {
                let (values, msg) = float_list_to_vec3(&v, Point3f::new);
                if let Some(msg) = msg {
                    warn!("point3d_param {}", msg);
                }
                Param::Point3d(name, values)
            },
        ))
    }

    /// Parse rule `vector3d_param`.
    fn vector3d_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [ident(name), float_list_expr(v)] => {
                let (values, msg) = float_list_to_vec3(&v, Vector3f::new);
                if let Some(msg) = msg {
                    warn!("vector3d_param {}", msg);
                }
                Param::Vector3d(name, values)
            },
        ))
    }

    /// Parse rule `normal3d_param`.
    fn normal3d_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [ident(name), float_list_expr(v)] => {
                let (values, msg) = float_list_to_vec3(&v, Normal3f::new);
                if let Some(msg) = msg {
                    warn!("normal3d_param {}", msg);
                }
                Param::Normal3d(name, values)
            },
        ))
    }

    /// Parse rule `point2d_param`.
    fn point2d_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [ident(name), float_list_expr(v)] => {
                let (values, msg) = float_list_to_vec2(&v, Point2f::new);
                if let Some(msg) = msg {
                    warn!("point2d_param {}", msg);
                }
                Param::Point2d(name, values)
            },
        ))
    }

    /// Parse rule `vector2d_param`.
    fn vector2d_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [ident(name), float_list_expr(v)] => {
                let (values, msg) = float_list_to_vec2(&v, Vector2f::new);
                if let Some(msg) = msg {
                    warn!("vector2d_param {}", msg);
                }
                Param::Vector2d(name, values)
            },
        ))
    }

    /// Parse rule `string_param`.
    fn string_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [ident(name), quoted_str_expr(v)] => Param::String(name, vec![v]),
            [ident(name), quoted_str_list_expr(v)] => Param::String(name, v),
        ))
    }

    /// Parse rule `bool_param`.
    fn bool_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [ident(name), quoted_bool_expr(v)] => Param::Bool(name, vec![v]),
            [ident(name), quoted_bool_list_expr(v)] => Param::Bool(name, v),
        ))
    }

    /// Parse rule `float_param`.
    fn float_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [ident(name), float_expr(v)] => Param::Float(name, vec![v]),
            [ident(name), float_list_expr(v)] => Param::Float(name, v),
        ))
    }

    /// Parse rule `int_param`.
    fn int_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [ident(name), int_expr(v)] => Param::Int(name, vec![v]),
            [ident(name), int_list_expr(v)] => Param::Int(name, v),
        ))
    }

    /// Parse rule `spectrum_param`.
    fn spectrum_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [ident(name), float_list_expr(v)] => Param::SpectrumFloat(name, v),
            [ident(name), quoted_str_expr(v)] => Param::SpectrumFile(name, v),
        ))
    }

    /// Parse rule `blackbody_param`.
    fn blackbody_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [ident(name), float_list_expr(v)] => Param::Blackbody(name, v),
        ))
    }

    /// Parse rule `texture_param`.
    fn texture_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [ident(name), quoted_str_expr(v)] => Param::Texture(name, vec![v]),
            [ident(name), quoted_str_list_expr(v)] => Param::Texture(name, v),
        ))
    }

    /// Parse rule `colour_param`.
    fn colour_param(input: Node) -> Result<Param> {
        Ok(match_nodes!(input.into_children();
            [colour_type(t), ident(name), float_list_expr(v)] =>
                Param::Colour(name, t, v),
        ))
    }

    /// Parse rule `colour_type`.
    fn colour_type(input: Node) -> Result<ColourType> {
        match input.as_str() {
            "colour" | "color" | "rgb" => Ok(ColourType::RGB),
            "xyz" => Ok(ColourType::XYZ),
            s => Err(input.error(ColourTypeError(s))),
        }
    }

    /// Parse rule `float_list_expr`.
    fn float_list_expr(input: Node) -> Result<Vec<Float>> {
        Ok(match_nodes!(input.into_children();
            [float_expr(values)..] => values.collect(),
        ))
    }

    /// Parse rule `float_expr`.
    fn float_expr(input: Node) -> Result<Float> {
        Ok(match_nodes!(input.into_children();
            [float(v)] => v,
        ))
    }

    /// Parse rule `float`.
    fn float(input: Node) -> Result<Float> {
        input.as_str().parse::<Float>().map_err(|e| input.error(e))
    }

    /// Parse rule `int_list_expr`.
    fn int_list_expr(input: Node) -> Result<Vec<Int>> {
        Ok(match_nodes!(input.into_children();
            [int_expr(values)..] => values.collect(),
        ))
    }

    /// Parse rule `int_expr`.
    fn int_expr(input: Node) -> Result<Int> {
        Ok(match_nodes!(input.into_children();
            [int(v)] => v,
        ))
    }

    /// Parse rule `int`.
    fn int(input: Node) -> Result<Int> {
        input.as_str().parse::<Int>().map_err(|e| input.error(e))
    }

    /// Parse rule `quoted_bool_list_expr`.
    fn quoted_bool_list_expr(input: Node) -> Result<Vec<bool>> {
        Ok(match_nodes!(input.into_children();
            [quoted_bool_expr(values)..] => values.collect(),
        ))
    }

    /// Parse rule `quoted_bool_expr`.
    fn quoted_bool_expr(input: Node) -> Result<bool> {
        Ok(match_nodes!(input.into_children();
            [quoted_bool(v)] => v,
        ))
    }

    /// Parse rule `quoted_bool`.
    fn quoted_bool(input: Node) -> Result<bool> {
        Ok(match_nodes!(input.into_children();
            [bool(v)] => v,
        ))
    }

    /// Parse rule `bool`.
    fn bool(input: Node) -> Result<bool> {
        input.as_str().parse::<bool>().map_err(|e| input.error(e))
    }

    /// Parse rule `quoted_str_list_expr`.
    fn quoted_str_list_expr(input: Node) -> Result<Vec<String>> {
        Ok(match_nodes!(input.into_children();
            [quoted_str_expr(values)..] => values.collect(),
        ))
    }

    /// Parse rule `quoted_str_expr`.
    fn quoted_str_expr(input: Node) -> Result<String> {
        Ok(match_nodes!(input.into_children();
            [quoted_str(v)] => v,
        ))
    }

    /// Parse rule `quoted_str`.
    fn quoted_str(input: Node) -> Result<String> {
        Ok(match_nodes!(input.into_children();
            [str(v)] => v,
        ))
    }

    /// Parse rule `str`.
    fn str(input: Node) -> Result<String> {
        Ok(input.as_str().to_owned())
    }

    /// Parse rule `quoted_ident_expr`.
    fn quoted_ident_expr(input: Node) -> Result<String> {
        Ok(match_nodes!(input.into_children();
            [quoted_ident(v)] => v,
        ))
    }

    /// Parse rule `quoted_ident`.
    fn quoted_ident(input: Node) -> Result<String> {
        Ok(match_nodes!(input.into_children();
            [ident(v)] => v,
        ))
    }

    /// Parse rule `ident_expr`.
    #[allow(unused)]
    fn ident_expr(input: Node) -> Result<String> {
        Ok(match_nodes!(input.into_children();
            [ident(v)] => v,
        ))
    }

    /// Parse rule `ident`.
    fn ident(input: Node) -> Result<String> {
        Ok(input.as_str().to_owned())
    }
}

/// Reads a PBRT file format and calls the API wrapper functions.
///
/// * `abs_path`        - The absolute path to scene file or an `Include` file.
/// * `unparsed_params` - Trailing parameters to append to `Include` file. Use empty string for main scene file.
///                       *NOTE:* If the `Include` file does not account for this, errors will point to lines from the
///                       original files as though they are in the include file.
/// * `api`             - The PBRT API interface.
pub fn parse(abs_path: &str, unparsed_params: &str, api: &mut Api) -> std::result::Result<(), String> {
    // Get path to scene file's folder for resolving relative paths to includes, images, etc.
    let scene_path = api.cwd.clone();

    // Load input file and append trailing parameters if any.
    let mut unparsed_file = file_to_string(&abs_path)?;
    if !unparsed_params.is_empty() {
        unparsed_file.push(' '); // Add whitespace to account for no leading whitespace in unparsed_params.
        unparsed_file.push_str(unparsed_params);
    }

    // Parse the input file into `Nodes`.
    match PbrtParser::parse_with_userdata(Rule::pbrt, &unparsed_file, scene_path)
        // There should be a single root node in the parsed tree.
        .and_then(|inputs| inputs.single())
        // Consume the `Node` recursively into the final value.
        .and_then(PbrtParser::pbrt)
        // Process the and call API.
        .and_then(|pbrt| {
            pbrt.process(api);
            Ok(())
        }) {
        Ok(()) => Ok(()),
        Err(e) => Err(format!("{:}", e)),
    }
}

/// Converts a named parameter to a vector of some type that stores 3 values.
///
/// * `v`    - Slice containing floating point values.
/// * `new`  - Function used to construct elements of resulting vector.
fn float_list_to_vec3<T, F>(v: &[Float], new: F) -> (Vec<T>, Option<String>)
where
    F: Fn(Float, Float, Float) -> T,
{
    let mut msg: Option<String> = None;

    let n = v.len();
    if n % 3 != 0 {
        msg = Some("length is not divisible by 3".to_owned());
    }
    let res = (0..n).step_by(3).map(|i| new(v[i], v[i + 1], v[i + 2])).collect();
    (res, msg)
}

/// Converts a named parameter to a vector of some type that stores 2 values.
///
/// * `v`    - Slice containing floating point values.
/// * `new`  - Function used to construct elements of resulting vector.
fn float_list_to_vec2<T, F>(v: &[Float], new: F) -> (Vec<T>, Option<String>)
where
    F: Fn(Float, Float) -> T,
{
    let mut msg: Option<String> = None;

    let n = v.len();
    if n % 2 != 0 {
        msg = Some("length is not divisible by 2".to_owned());
    }
    let res = (0..n).step_by(2).map(|i| new(v[i], v[i + 1])).collect();
    (res, msg)
}
