//! Material Instance

#![allow(dead_code)]
use crate::core::material::*;
use crate::core::paramset::*;
use std::sync::Arc;

/// Stores a reference to a named material and its parameters.
#[derive(Clone)]
pub struct MaterialInstance {
    /// The name.
    pub name: String,

    /// Reference to the material.
    pub material: ArcMaterial,

    /// The parameters.
    pub params: Arc<ParamSet>,
}

impl MaterialInstance {
    /// Create a new `MaterialInstance`.
    ///
    /// * `name`     - The name.
    /// * `material` - Reference to the material.
    /// * `params`   - Parameters.
    pub fn new(name: &str, material: ArcMaterial, params: Arc<ParamSet>) -> Self {
        Self {
            name: String::from(name),
            material: material.clone(),
            params: params.clone(),
        }
    }
}
