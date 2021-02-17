//! PBRT File Parser

#![allow(dead_code)]

use crate::core::api::*;
use crate::core::fileutil::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;
use pest::iterators::*;
use pest::Parser;
use std::fs;
use std::result::Result;

/// The `pest` parser generated from a grammar.
#[derive(Parser)]
#[grammar = "core/parsers/pbrt/grammar.pest"]
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
        match next_pair.as_rule() {
            Rule::empty_stmt => (),   // Ignore
            Rule::comment_stmt => (), // Ignore
            Rule::include_stmt => {
                let mut inner_rules = next_pair.into_inner();
                self.parse_include_stmt(&mut inner_rules, api);
            }
            Rule::option_stmt => {
                let mut inner_rules = next_pair.into_inner();
                self.parse_option_stmt(&mut inner_rules, api);
            }
            Rule::scene_stmt => (),
            Rule::block_stmt => {
                let mut inner_rules = next_pair.into_inner();
                self.parse_block_stmt(&mut inner_rules, api);
            }
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
    fn parse_option_stmt(&self, pairs: &mut Pairs<Rule>, _api: &mut Api) {
        let next_pair = pairs.next().unwrap();
        let rule = next_pair.as_rule();
        let mut inner_rules = next_pair.into_inner();

        let name = self.parse_quoted_str(&mut inner_rules);

        match rule {
            Rule::accelerator_stmt => (),
            Rule::camera_stmt => (),
            Rule::film_stmt => {
                let param_list = inner_rules.next().unwrap();
                let params = self.parse_param_list(param_list.into_inner());
                debug!("Film '{}' {:}", name, params);
                // TODO create film.
            }
            Rule::filter_stmt => (),
            Rule::integrator_stmt => (),
            Rule::make_named_medium_stmt => (),
            Rule::sampler_stmt => (),
            _ => unreachable!(),
        }
    }

    /// Parse a `scene_stmt` rule of the grammar and call the API.
    ///
    /// * `pairs` - The inner token pairs for matched `scene_stmt` rule.
    /// * `api`   - The PBRT API interface.
    fn parse_scene_stmt(&self, pairs: &mut Pairs<Rule>, _api: &mut Api) {
        let next_pair = pairs.next().unwrap();
        match next_pair.as_rule() {
            Rule::area_light_source_stmt => (),
            Rule::light_source_stmt => (),
            Rule::make_named_material_stmt => (),
            Rule::material_stmt => (),
            Rule::shape_stmt => (),
            Rule::texture_stmt => (),
            Rule::named_material_stmt => (),
            Rule::object_instance_stmt => (),
            Rule::reverse_orientation_stmt => (),
            Rule::medium_interface_stmt => (),
            Rule::ctm_stmt => (),
            Rule::active_transform_stmt => (),
            Rule::transform_type => (),
            _ => unreachable!(),
        }
    }

    /// Parse a `ctm_stmt` rule of the grammar and call the API.
    ///
    /// * `pairs` - The inner token pairs for matched `ctm_stmt` rule.
    /// * `api`   - The PBRT API interface.
    fn parse_ctm_stmt(&self, pairs: &mut Pairs<Rule>, _api: &mut Api) {
        let next_pair = pairs.next().unwrap();
        match next_pair.as_rule() {
            Rule::identity_stmt => (),
            Rule::translate_stmt => (),
            Rule::scale_stmt => (),
            Rule::rotate_stmt => (),
            Rule::look_at_stmt => (),
            Rule::coordinate_system_stmt => (),
            Rule::coord_sys_transform_stmt => (),
            Rule::transform_stmt => (),
            Rule::concat_transform_stmt => (),
            Rule::transform_times_stmt => (),
            _ => unreachable!(),
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
            Rule::point3d_param => (),
            Rule::vector3d_param => (),
            Rule::normal3d_param => (),
            Rule::point2d_param => (),
            Rule::vector2d_param => (),
            Rule::string_param => (),
            Rule::bool_param => (),
            Rule::float_param => (),
            Rule::int_param => self.parse_int_param(&mut inner_rules, params),
            Rule::colour_param => (),
            Rule::spectrum_param => (),
            Rule::blackbody_param => (),
            _ => unreachable!(),
        }
    }

    /// Parse an `int_param` rule of the grammar and add parameter to a
    /// `ParamSet`.
    ///
    /// * `pairs`  - The inner token pairs for matched `int_list` rule.
    /// * `params` - The `ParamSet` to update.
    fn parse_int_param(&self, pairs: &mut Pairs<Rule>, params: &mut ParamSet) {
        let int_type = pairs.next().unwrap().as_str();
        assert!(int_type == "integer");

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

    /// Parse an `int_expr` or `int` rule of the grammar and return an `Int`.
    ///
    /// * `pairs`  - The inner token pairs for matched `int_expr` or `int` rule.
    fn parse_int(&self, pair: Pair<Rule>) -> Int {
        // Parse string to int. The unwrap shouldn't fail if our pest
        // grammar is correct.
        pair.as_str().parse::<Int>().unwrap()
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
