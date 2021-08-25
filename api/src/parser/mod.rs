//! PBRT File Parser

#![allow(dead_code)]

use super::*;
use core::fileutil::*;
use pest::iterators::*;
use pest::Parser;
use std::fs;
use std::result::Result;

/// The `pest` parser generated from a grammar.
#[derive(Parser)]
#[grammar = "parser/grammar.pest"]
struct PbrtParser;

/// PBRT File Format Parser.
pub struct PbrtFileParser {
    /// Path to the file to parse.
    file_path: String,

    /// Parent path for navigating to includes.
    parent_path: String,
}

impl PbrtFileParser {
    /// Returns a new instance of `PbrtFileParser`.
    ///
    /// * `path` - File path.
    pub fn new(path: &str) -> Self {
        if let Some(parent) = parent_path(path) {
            Self {
                file_path: String::from(path),
                parent_path: parent,
            }
        } else {
            // We were passed the root path itself which is not a file.
            panic!("Invalid path '{}'", path);
        }
    }

    /// Reads a PBRT file format and calls the API wrapper functions.
    ///
    /// * `api`  - The PBRT API interface.
    pub fn parse(&self, api: &mut Api) -> Result<(), String> {
        // Load the file and parse the `file` rule.
        let unparsed_file = file_to_string(&self.file_path)?;
        let pbrt = self.parse_pbrt_rule(&unparsed_file)?;

        // Parse all the `stmt` rules.
        for pair in pbrt.into_inner() {
            match pair.as_rule() {
                Rule::stmt => {
                    let mut inner_rules = pair.into_inner();
                    self.parse_stmt_rule(&mut inner_rules, api);
                }
                Rule::EOI => (), // Done
                _ => unreachable!(),
            }
        }

        Ok(())
    }

    /// Parse the initial `pbrt` rule of the grammar and return the resulting token
    /// pairs for remaining rules.
    ///
    /// * `unparsed_file` - Contents of the file to parse.
    fn parse_pbrt_rule<'a>(&self, unparsed_file: &'a str) -> Result<Pair<'a, Rule>, String> {
        match PbrtParser::parse(Rule::pbrt, &unparsed_file) {
            Ok(mut pairs) => Ok(pairs.next().unwrap()),
            Err(err) => Err(format!("Error parsing pbrt rule. {}", err)),
        }
    }

    /// Parse a `stmt` rule of the grammar and call the API.
    ///
    /// * `pairs` - The inner token pairs for matched `stmt` rule.
    fn parse_stmt_rule(&self, pairs: &mut Pairs<Rule>, api: &mut Api) {
        let next_pair = pairs.next().unwrap();
        let rule = next_pair.as_rule();
        let mut inner_rules = next_pair.into_inner();

        match rule {
            Rule::empty_stmt => (),   // Ignore
            Rule::comment_stmt => (), // Ignore
            Rule::include_stmt => self.parse_include_stmt(&mut inner_rules, api),
            Rule::option_stmt => self.parse_option_stmt(&mut inner_rules, api),
            Rule::scene_stmt => self.parse_scene_stmt(&mut inner_rules, api),
            Rule::block_stmt => self.parse_block_stmt(&mut inner_rules, api),
            Rule::ctm_stmt => self.parse_ctm_stmt(&mut inner_rules, api),
            _ => unreachable!(),
        }
    }

    /// Parse a `include_stmt` rule of the grammar and call the API.
    /// This will create a new parser to parse the included file and
    /// parse it entirely while calling the API before returning.
    ///
    /// * `pairs` - The inner token pairs for matched `include_stmt` rule.
    /// * `api`   - The PBRT API interface.
    fn parse_include_stmt(&self, pairs: &mut Pairs<Rule>, api: &mut Api) {
        let next_pair = pairs.next().unwrap();
        match next_pair.as_rule() {
            Rule::quoted_str_expr => {
                let mut inner_rules = next_pair.into_inner();
                let mut path = self.parse_quoted_str(&mut inner_rules);
                debug!("Include: '{}'", path);

                if is_relative_path(&path) {
                    // Path is relative to the parent path of the file being parsed.
                    path = self.parent_path.clone() + "/" + &path;
                }

                let parser = Self::new(&path);
                match parser.parse(api) {
                    Ok(()) => debug!("Finished parsing include '{}'", path),
                    Err(err) => error!("{}", err),
                }
            }
            _ => unreachable!(),
        }
    }

    /// Parse a `block_stmt` rule of the grammar and call the API.
    ///
    /// * `pairs` - The inner token pairs for matched `block_stmt` rule.
    /// * `api`   - The PBRT API interface.
    fn parse_block_stmt(&self, pairs: &mut Pairs<Rule>, api: &mut Api) {
        let next_pair = pairs.next().unwrap();
        match next_pair.as_rule() {
            Rule::world_begin_stmt => api.pbrt_world_begin(),
            Rule::world_end_stmt => api.pbrt_world_end(),
            Rule::attribute_begin_stmt => api.pbrt_attribute_begin(),
            Rule::attribute_end_stmt => api.pbrt_attribute_end(),
            Rule::object_begin_stmt => {
                let mut inner_rules = next_pair.into_inner();
                let str = self.parse_quoted_str(&mut inner_rules);
                debug!("ObjectBegin: '{}'", str);
                api.pbrt_object_begin(str);
            }
            Rule::object_end_stmt => api.pbrt_object_end(),
            _ => unreachable!(),
        }
    }

    /// Parse an `option_stmt` rule of the grammar and call the API.
    ///
    /// * `pairs` - The inner token pairs for matched `option_stmt` rule.
    /// * `api`   - The PBRT API interface.
    fn parse_option_stmt(&self, pairs: &mut Pairs<Rule>, api: &mut Api) {
        let next_pair = pairs.next().unwrap();
        let rule = next_pair.as_rule();
        let mut inner_rules = next_pair.into_inner();

        // TODO Create the appropriate objects via the API.
        match rule {
            Rule::accelerator_stmt => {
                self.parse_named_param_list(&mut inner_rules, "Accelerator", api)
            }
            Rule::camera_stmt => self.parse_named_param_list(&mut inner_rules, "Camera", api),
            Rule::film_stmt => self.parse_named_param_list(&mut inner_rules, "Film", api),
            Rule::filter_stmt => self.parse_named_param_list(&mut inner_rules, "Filter", api),
            Rule::integrator_stmt => {
                self.parse_named_param_list(&mut inner_rules, "Integrator", api)
            }
            Rule::make_named_medium_stmt => {
                self.parse_named_param_list(&mut inner_rules, "MakeNamedMedium", api)
            }
            Rule::sampler_stmt => self.parse_named_param_list(&mut inner_rules, "Sampler", api),
            _ => unreachable!(),
        }
    }

    /// Parse a `scene_stmt` rule of the grammar and call the API.
    ///
    /// * `pairs` - The inner token pairs for matched `scene_stmt` rule.
    /// * `api`   - The PBRT API interface.
    fn parse_scene_stmt(&self, pairs: &mut Pairs<Rule>, api: &mut Api) {
        let next_pair = pairs.next().unwrap();
        let rule = next_pair.as_rule();
        let mut inner_rules = next_pair.into_inner();

        // TODO Create the appropriate objects via the API.
        match rule {
            Rule::area_light_source_stmt => {
                self.parse_named_param_list(&mut inner_rules, "AreaLightSource", api)
            }
            Rule::light_source_stmt => {
                self.parse_named_param_list(&mut inner_rules, "LightSource", api)
            }
            Rule::make_named_material_stmt => {
                self.parse_named_param_list(&mut inner_rules, "MakeNamedMaterial", api)
            }
            Rule::material_stmt => self.parse_named_param_list(&mut inner_rules, "Material", api),
            Rule::shape_stmt => self.parse_named_param_list(&mut inner_rules, "Shape", api),
            Rule::texture_stmt => {
                let name = self.parse_quoted_str(&mut inner_rules);
                let texture_type = self.parse_quoted_str(&mut inner_rules);
                let texture_name = self.parse_quoted_str(&mut inner_rules);
                let params = inner_rules.next().map_or(ParamSet::new(), |param_list| {
                    self.parse_param_list(param_list.into_inner())
                });
                debug!(
                    "Texture: '{}', '{}', '{}' {:}",
                    name, texture_type, texture_name, params
                );
                api.pbrt_texture(name, texture_type, texture_name, &params)
            }
            Rule::named_material_stmt => {
                self.parse_named_stmt(&mut inner_rules, api, "NamedMaterial")
            }
            Rule::object_instance_stmt => {
                self.parse_named_stmt(&mut inner_rules, api, "ObjectInstance")
            }
            Rule::reverse_orientation_stmt => api.pbrt_reverse_orientation(),
            Rule::medium_interface_stmt => {
                let inside_medium = inner_rules.next().unwrap().as_str().to_string();
                let outside_medium = inner_rules.next().unwrap().as_str().to_string();
                debug!("MediumInterface: '{}', '{}'", inside_medium, outside_medium);
                api.pbrt_medium_interface(inside_medium, outside_medium);
            }
            Rule::ctm_stmt => self.parse_ctm_stmt(&mut inner_rules, api),
            Rule::active_transform_stmt => {
                let time = inner_rules.next().unwrap().as_str();
                debug!("ActiveTransform: '{}'", time);
                match time {
                    "StartTime" => api.pbrt_active_transform_start_time(),
                    "EndTime" => api.pbrt_active_transform_end_time(),
                    "All" => api.pbrt_active_transform_all(),
                    _ => warn!("Ignoring invalid ActiveTransform time '{}'", time),
                }
            }
            _ => unreachable!(),
        }
    }

    /// Parse a `named_material_stmt`/`object_instance_stmt` rule of the grammar and call the API.
    ///
    /// * `pairs`     - The inner token pairs for matched rule.
    /// * `api`       - The PBRT API interface.
    /// * `stmt_type` - 'NamedMaterial | ObjectInstance'
    fn parse_named_stmt(&self, pairs: &mut Pairs<Rule>, api: &mut Api, stmt_type: &str) {
        let next_pair = pairs.next().unwrap();
        match next_pair.as_rule() {
            Rule::quoted_ident_expr => {
                let mut inner_rules = next_pair.into_inner();
                let name = self.parse_quoted_ident(&mut inner_rules);
                debug!("{}: '{}'", stmt_type, name);
                match stmt_type {
                    "NamedMaterial" => api.pbrt_named_material(name),
                    "ObjectInstance" => api.pbrt_object_instance(name),
                    _ => unreachable!(),
                }
            }
            _ => unreachable!(),
        }
    }

    /// Parse a `ctm_stmt` rule of the grammar and call the API.
    ///
    /// * `pairs` - The inner token pairs for matched `ctm_stmt` rule.
    /// * `api`   - The PBRT API interface.
    fn parse_ctm_stmt(&self, pairs: &mut Pairs<Rule>, api: &mut Api) {
        let next_pair = pairs.next().unwrap();
        match next_pair.as_rule() {
            Rule::identity_stmt => {
                debug!("Identity");
                api.pbrt_identity();
            }
            Rule::translate_stmt => {
                let mut inner_rules = next_pair.into_inner();
                let x = self.parse_float(inner_rules.next().unwrap());
                let y = self.parse_float(inner_rules.next().unwrap());
                let z = self.parse_float(inner_rules.next().unwrap());
                debug!("Translate: [{}, {}, {}]", x, y, z);
                api.pbrt_translate(x, y, z);
            }
            Rule::scale_stmt => {
                let mut inner_rules = next_pair.into_inner();
                let x = self.parse_float(inner_rules.next().unwrap());
                let y = self.parse_float(inner_rules.next().unwrap());
                let z = self.parse_float(inner_rules.next().unwrap());
                debug!("Scale: [{}, {}, {}]", x, y, z);
                api.pbrt_scale(x, y, z);
            }
            Rule::rotate_stmt => {
                let mut inner_rules = next_pair.into_inner();
                let angle = self.parse_float(inner_rules.next().unwrap());
                let x = self.parse_float(inner_rules.next().unwrap());
                let y = self.parse_float(inner_rules.next().unwrap());
                let z = self.parse_float(inner_rules.next().unwrap());
                debug!("Rotate: {}, [{}, {}, {}]", angle, x, y, z);
                api.pbrt_rotate(angle, x, y, z);
            }
            Rule::look_at_stmt => {
                let mut inner_rules = next_pair.into_inner();
                let ex = self.parse_float(inner_rules.next().unwrap());
                let ey = self.parse_float(inner_rules.next().unwrap());
                let ez = self.parse_float(inner_rules.next().unwrap());
                let lx = self.parse_float(inner_rules.next().unwrap());
                let ly = self.parse_float(inner_rules.next().unwrap());
                let lz = self.parse_float(inner_rules.next().unwrap());
                let ux = self.parse_float(inner_rules.next().unwrap());
                let uy = self.parse_float(inner_rules.next().unwrap());
                let uz = self.parse_float(inner_rules.next().unwrap());
                debug!(
                    "LookAt: [{}, {}, {}], [{}, {}, {}], [{}, {}, {}]",
                    ex, ey, ez, lx, ly, lz, ux, uy, uz
                );
                api.pbrt_look_at(ex, ey, ez, lx, ly, lz, ux, uy, uz);
            }
            Rule::coordinate_system_stmt => {
                let mut inner_rules = next_pair.into_inner();
                let name = self.parse_quoted_ident(&mut inner_rules);
                debug!("CoordinateSystem: '{}'", name);
                api.pbrt_coordinate_system(name);
            }
            Rule::coord_sys_transform_stmt => {
                let mut inner_rules = next_pair.into_inner();
                let name = self.parse_quoted_ident(&mut inner_rules);
                debug!("CoordSysTransform: '{}'", name);
                api.pbrt_coord_sys_transform(name);
            }
            Rule::transform_stmt => {
                let tr = self.parse_float_list(next_pair.into_inner());
                assert!(
                    tr.len() == 16,
                    "float_list in transform_stmt not of len 16."
                );
                debug!("Transform: {:?}", tr);
                api.pbrt_transform(&[
                    tr[0], tr[1], tr[2], tr[3], tr[4], tr[5], tr[6], tr[7], tr[8], tr[9], tr[10],
                    tr[11], tr[11], tr[12], tr[13], tr[14],
                ]);
            }
            Rule::concat_transform_stmt => {
                let tr = self.parse_float_list(next_pair.into_inner());
                assert!(
                    tr.len() == 16,
                    "float_list in concat_transform_stmt not of len 16."
                );
                debug!("ConcatTransform: {:?}", tr);
                api.pbrt_concat_transform(&[
                    tr[0], tr[1], tr[2], tr[3], tr[4], tr[5], tr[6], tr[7], tr[8], tr[9], tr[10],
                    tr[11], tr[11], tr[12], tr[13], tr[14],
                ]);
            }
            Rule::transform_times_stmt => {
                let mut inner_rules = next_pair.into_inner();
                let start = self.parse_float(inner_rules.next().unwrap());
                let end = self.parse_float(inner_rules.next().unwrap());
                debug!("TransformTimes: {}, {}", start, end);
                api.pbrt_transform_times(start, end);
            }
            _ => unreachable!(),
        }
    }

    /// Parse a named `param_list`.
    ///
    /// * `pairs`       - The inner token pairs for matched `param_list` rule.
    /// * `option_name` - The name of the option.
    /// * `api`         - The PBRT API interface.
    fn parse_named_param_list(&self, pairs: &mut Pairs<Rule>, option_name: &str, api: &mut Api) {
        let name = self.parse_quoted_str(pairs);
        let params = pairs.next().map_or(ParamSet::new(), |param_list| {
            self.parse_param_list(param_list.into_inner())
        });

        debug!("{} '{}' {:}", option_name, name, params);
        match option_name {
            "Accelerator" => api.pbrt_accelerator(name, &params),
            "Camera" => api.pbrt_camera(name, &params),
            "Film" => api.pbrt_film(name, &params),
            "Filter" => api.pbrt_pixel_filter(name, &params),
            "Integrator" => api.pbrt_integrator(name, &params),
            "MakeNamedMedium" => api.pbrt_make_named_medium(name, &params),
            "Sampler" => api.pbrt_sampler(name, &params),
            "AreaLightSource" => api.pbrt_area_light_source(name, &params),
            "LightSource" => api.pbrt_light_source(name, &params),
            "MakeNamedMaterial" => api.pbrt_make_named_material(name, &params),
            "Material" => api.pbrt_material(name, &params),
            "Shape" => api.pbrt_shape(name, &params),
            _ => warn!("'{}' not supported", option_name),
        }
    }

    /// Parse a `param_list` rule of the grammar and return a `ParamSet`.
    ///
    /// * `pairs` - The inner token pairs for matched `param_list` rule.
    fn parse_param_list(&self, pairs: Pairs<Rule>) -> ParamSet {
        let mut paramset = ParamSet::new();

        for pair in pairs {
            let rule = pair.as_rule();
            let mut inner_rules = pair.into_inner();
            match rule {
                Rule::param => self.parse_param(&mut inner_rules, &mut paramset),
                Rule::comment => (), // Ignore
                _ => unreachable!(),
            }
        }

        paramset
    }

    /// Parse a `param` rule of the grammar and add parameter to a `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `param_list` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let next_pair = pairs.next().unwrap();
        let rule = next_pair.as_rule();
        let mut inner_rules = next_pair.into_inner();
        match rule {
            Rule::point3d_param => self.parse_point3d_param(&mut inner_rules, params),
            Rule::vector3d_param => self.parse_vector3d_param(&mut inner_rules, params),
            Rule::normal3d_param => self.parse_normal3d_param(&mut inner_rules, params),
            Rule::point2d_param => self.parse_point2d_param(&mut inner_rules, params),
            Rule::vector2d_param => self.parse_vector2d_param(&mut inner_rules, params),
            Rule::string_param => self.parse_string_param(&mut inner_rules, params),
            Rule::bool_param => self.parse_bool_param(&mut inner_rules, params),
            Rule::float_param => self.parse_float_param(&mut inner_rules, params),
            Rule::int_param => self.parse_int_param(&mut inner_rules, params),
            Rule::colour_param => self.parse_colour_param(&mut inner_rules, params),
            Rule::spectrum_param => self.parse_spectrum_param(&mut inner_rules, params),
            Rule::blackbody_param => self.parse_blackbody_param(&mut inner_rules, params),
            Rule::texture_param => self.parse_texture_param(&mut inner_rules, params),
            _ => unreachable!(),
        }
    }

    /// Parse an `point3d_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `point3d_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_point3d_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "point3" || param_type == "point");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        let list = match value.as_rule() {
            Rule::float_list_expr => self.parse_float_list(value.into_inner()),
            _ => unreachable!(),
        };

        let n = list.len();
        if n % 3 != 0 {
            warn!("point3d_param '{}' length is not divisible by 3", ident);
        }

        let values: Vec<Point3f> = (0..n)
            .step_by(3)
            .map(|i| Point3f::new(list[i], list[i + 1], list[i + 2]))
            .collect();
        params.add_point3f(ident, &values);
    }

    /// Parse an `vector3d_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `vector3d_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_vector3d_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "vector3" || param_type == "vector");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        let list = match value.as_rule() {
            Rule::float_list_expr => self.parse_float_list(value.into_inner()),
            _ => unreachable!(),
        };

        let n = list.len();
        if n % 3 != 0 {
            warn!("vector3d_param '{}' length is not divisible by 3", ident);
        }

        let values: Vec<Vector3f> = (0..n)
            .step_by(3)
            .map(|i| Vector3f::new(list[i], list[i + 1], list[i + 2]))
            .collect();
        params.add_vector3f(ident, &values);
    }

    /// Parse an `normal3d_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `normal3d_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_normal3d_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "normal3" || param_type == "normal");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        let list = match value.as_rule() {
            Rule::float_list_expr => self.parse_float_list(value.into_inner()),
            _ => unreachable!(),
        };

        let n = list.len();
        if n % 3 != 0 {
            warn!("normal3d_param '{}' length is not divisible by 3", ident);
        }

        let values: Vec<Normal3f> = (0..n)
            .step_by(3)
            .map(|i| Normal3f::new(list[i], list[i + 1], list[i + 2]))
            .collect();
        params.add_normal3f(ident, &values);
    }

    /// Parse an `point2d_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `point2d_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_point2d_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "point2");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        let list = match value.as_rule() {
            Rule::float_list_expr => self.parse_float_list(value.into_inner()),
            _ => unreachable!(),
        };

        let n = list.len();
        if n % 2 != 0 {
            warn!("point2d_param '{}' length is not divisible by 3", ident);
        }

        let values: Vec<Point2f> = (0..n)
            .step_by(2)
            .map(|i| Point2f::new(list[i], list[i + 1]))
            .collect();
        params.add_point2f(ident, &values);
    }

    /// Parse an `vector2d_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `vector2d_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_vector2d_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "vector2");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        let list = match value.as_rule() {
            Rule::float_list_expr => self.parse_float_list(value.into_inner()),
            _ => unreachable!(),
        };

        let n = list.len();
        if n % 2 != 0 {
            warn!("vector2d_param '{}' length is not divisible by 3", ident);
        }

        let values: Vec<Vector2f> = (0..n)
            .step_by(2)
            .map(|i| Vector2f::new(list[i], list[i + 1]))
            .collect();
        params.add_vector2f(ident, &values);
    }

    /// Parse an `string_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `string_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_string_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "string");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        let list: Vec<String> = match value.as_rule() {
            Rule::quoted_str_expr => {
                let mut inner_rules = value.into_inner();
                vec![self.parse_quoted_str(&mut inner_rules)]
            }
            Rule::quoted_str_list_expr => self.parse_quoted_str_list(value.into_inner()),
            _ => unreachable!(),
        };
        params.add_string(ident, &list);
    }

    /// Parse an `bool_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `bool_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_bool_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "bool");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        let list: Vec<bool> = match value.as_rule() {
            Rule::quoted_bool_expr => {
                let mut inner_rules = value.into_inner();
                vec![self.parse_quoted_bool(&mut inner_rules)]
            }
            Rule::quoted_bool_list_expr => self.parse_quoted_bool_list(value.into_inner()),
            _ => unreachable!(),
        };
        params.add_bool(ident, &list);
    }

    /// Parse an `float_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `float_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_float_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "float");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        let list = match value.as_rule() {
            Rule::float_expr => vec![self.parse_float(value)],
            Rule::float_list_expr => self.parse_float_list(value.into_inner()),
            _ => unreachable!(),
        };
        params.add_float(ident, &list);
    }

    /// Parse an `int_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `int_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_int_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "integer");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        let list = match value.as_rule() {
            Rule::int_expr => vec![self.parse_int(value)],
            Rule::int_list_expr => self.parse_int_list(value.into_inner()),
            _ => unreachable!(),
        };
        params.add_int(ident, &list);
    }

    /// Parse an `colour_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `colour_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_colour_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "color" || param_type == "rgb" || param_type == "xyz");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        let list = match value.as_rule() {
            Rule::float_list_expr => self.parse_float_list(value.into_inner()),
            _ => unreachable!(),
        };

        // ParamSet does additional validation.
        if param_type == "color" || param_type == "rgb" {
            params.add_rgb_spectrum(ident, &list);
        } else {
            params.add_xyz_spectrum(ident, &list);
        }
    }

    /// Parse an `spectrum_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `spectrum_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_spectrum_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "spectrum");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        match value.as_rule() {
            Rule::float_list_expr => {
                // ParamSet does additional validation.
                let list = self.parse_float_list(value.into_inner());
                params.add_sampled_spectrum(ident, &list);
            }
            Rule::quoted_str_expr => {
                let mut inner_rules = value.into_inner();
                let filename = self.parse_quoted_str(&mut inner_rules);
                params.add_sampled_spectrum_files(ident, &[filename]);
            }
            _ => unreachable!(),
        };
    }

    /// Parse an `blackbody_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `blackbody_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_blackbody_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "blackbody");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        let list = match value.as_rule() {
            Rule::float_list_expr => self.parse_float_list(value.into_inner()),
            _ => unreachable!(),
        };

        // ParamSet does additional validation.
        params.add_blackbody_spectrum(ident, &list);
    }

    /// Parse an `texture_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `texture_param` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_texture_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let param_type = pairs.next().unwrap().as_str();
        assert!(param_type == "texture");

        let ident = pairs.next().unwrap().as_str();
        let value = pairs.next().unwrap();
        assert!(pairs.next().is_none());

        let list: Vec<String> = match value.as_rule() {
            Rule::quoted_str_expr => {
                let mut inner_rules = value.into_inner();
                vec![self.parse_quoted_str(&mut inner_rules)]
            }
            Rule::quoted_str_list_expr => self.parse_quoted_str_list(value.into_inner()),
            _ => unreachable!(),
        };
        params.add_texture(ident, &list);
    }

    /// Parse an `float_list_expr` rule of the grammar and return a `Vec<Float>`.
    ///
    /// * `pairs`  - The inner token pairs for matched `float_list_expr` rule.
    fn parse_float_list(&self, pairs: Pairs<Rule>) -> Vec<Float> {
        let mut v: Vec<Float> = vec![];
        for pair in pairs {
            match pair.as_rule() {
                Rule::float => v.push(self.parse_float(pair)),
                _ => unreachable!(),
            }
        }
        v
    }

    /// Parse an `int_list_expr` rule of the grammar and return a `Vec<Int>`.
    ///
    /// * `pairs`  - The inner token pairs for matched `int_list_expr` rule.
    fn parse_int_list(&self, pairs: Pairs<Rule>) -> Vec<Int> {
        let mut v: Vec<Int> = vec![];
        for pair in pairs {
            match pair.as_rule() {
                Rule::int => v.push(self.parse_int(pair)),
                _ => unreachable!(),
            }
        }
        v
    }

    /// Parse an `quoted_str_list_expr` rule of the grammar and return a `Vec<String>`.
    ///
    /// * `pairs`  - The inner token pairs for matched `quoted_str_list_expr` rule.
    fn parse_quoted_str_list(&self, pairs: Pairs<Rule>) -> Vec<String> {
        let mut v: Vec<String> = vec![];
        for pair in pairs {
            match pair.as_rule() {
                Rule::quoted_str => {
                    let mut inner_rules = pair.into_inner();
                    v.push(self.parse_str(&mut inner_rules));
                }
                _ => unreachable!(),
            }
        }
        v
    }

    /// Parse an `quoted_bool_list_expr` rule of the grammar and return a `Vec<String>`.
    ///
    /// * `pairs`  - The inner token pairs for matched `quoted_bool_list_expr` rule.
    fn parse_quoted_bool_list(&self, pairs: Pairs<Rule>) -> Vec<bool> {
        let mut v: Vec<bool> = vec![];
        for pair in pairs {
            match pair.as_rule() {
                Rule::quoted_bool => {
                    let mut inner_rules = pair.into_inner();
                    v.push(self.parse_bool(&mut inner_rules));
                }
                _ => unreachable!(),
            }
        }
        v
    }

    /// Parse a `quoted_str` rule of the grammar and return the unquoted
    /// `String` value.
    ///
    /// * `pairs`  - The inner token pairs for matched `quoted_str` rule.
    fn parse_quoted_str(&self, pairs: &mut Pairs<Rule>) -> String {
        let next_pair = pairs.next().unwrap();
        match next_pair.as_rule() {
            Rule::quoted_str => {
                let mut inner_rules = next_pair.into_inner();
                self.parse_str(&mut inner_rules)
            }
            _ => unreachable!(),
        }
    }

    /// Parse a `quoted_ident` rule of the grammar and return the unquoted
    /// `String` value.
    ///
    /// * `pairs`  - The inner token pairs for matched `quoted_ident` rule.
    fn parse_quoted_ident(&self, pairs: &mut Pairs<Rule>) -> String {
        let next_pair = pairs.next().unwrap();
        match next_pair.as_rule() {
            Rule::quoted_ident => {
                let mut inner_rules = next_pair.into_inner();
                self.parse_ident(&mut inner_rules)
            }
            _ => unreachable!(),
        }
    }

    /// Parse a `quoted_bool` rule of the grammar and return the unquoted
    /// `bool` value.
    ///
    /// * `pairs`  - The inner token pairs for matched `quoted_bool` rule.
    fn parse_quoted_bool(&self, pairs: &mut Pairs<Rule>) -> bool {
        let next_pair = pairs.next().unwrap();
        match next_pair.as_rule() {
            Rule::quoted_bool => {
                let mut inner_rules = next_pair.into_inner();
                self.parse_bool(&mut inner_rules)
            }
            _ => unreachable!(),
        }
    }

    /// Parse an `float_expr` or `float` rule of the grammar and return an `Float`.
    ///
    /// * `pairs`  - The inner token pairs for matched `float_expr` or `float` rule.
    fn parse_float(&self, pair: Pair<Rule>) -> Float {
        // Parse string to float. The unwrap shouldn't fail if our pest
        // grammar is correct.
        let s = match pair.as_rule() {
            Rule::float_expr => {
                let mut inner_rules = pair.into_inner();
                inner_rules.next().unwrap().as_str()
            }
            Rule::float => pair.as_str(),
            _ => unreachable!(),
        };
        s.parse::<Float>().unwrap()
    }

    /// Parse an `int_expr` or `int` rule of the grammar and return an `Int`.
    ///
    /// * `pairs`  - The inner token pairs for matched `int_expr` or `int` rule.
    fn parse_int(&self, pair: Pair<Rule>) -> Int {
        // Parse string to int. The unwrap shouldn't fail if our pest
        // grammar is correct.
        let s = match pair.as_rule() {
            Rule::int_expr => {
                let mut inner_rules = pair.into_inner();
                inner_rules.next().unwrap().as_str()
            }
            Rule::int => pair.as_str(),
            _ => unreachable!(),
        };
        s.parse::<Int>().unwrap()
    }

    /// Parse a `str` rule of the grammar and return the `String` value.
    ///
    /// * `pairs`  - The inner token pairs for matched `str` rule.
    fn parse_str(&self, pairs: &mut Pairs<Rule>) -> String {
        let next_pair = pairs.next().unwrap();
        match next_pair.as_rule() {
            Rule::str => String::from(next_pair.as_str()),
            _ => unreachable!(),
        }
    }

    /// Parse an `ident` rule of the grammar and return the `String` value.
    ///
    /// * `pairs`  - The inner token pairs for matched `ident` rule.
    fn parse_ident(&self, pairs: &mut Pairs<Rule>) -> String {
        let next_pair = pairs.next().unwrap();
        match next_pair.as_rule() {
            Rule::ident => String::from(next_pair.as_str()),
            _ => unreachable!(),
        }
    }

    /// Parse a `bool` rule of the grammar and return a `bool`.
    ///
    /// * `pairs`  - The inner token pairs for matched `bool` rule.
    fn parse_bool(&self, pairs: &mut Pairs<Rule>) -> bool {
        let next_pair = pairs.next().unwrap();
        match next_pair.as_rule() {
            Rule::bool => next_pair.as_str().parse::<bool>().unwrap(),
            _ => unreachable!(),
        }
    }
}

/// Read the entire file and return its contents as a String.
///
/// * `path` - Path to file.
fn file_to_string(path: &str) -> Result<String, String> {
    match fs::read_to_string(path) {
        Ok(s) => Ok(s),
        _ => Err(format!("Error reading file '{}'", path)),
    }
}
