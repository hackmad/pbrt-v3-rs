//! Ray intersection acceleration data structures.

#![allow(dead_code)]
mod bvh;
mod kd_tree;

// Re-export
pub use bvh::*;
pub use kd_tree::*;

use crate::core::paramset::ParamSet;
use crate::core::primitive::ArcPrimitive;

/// Stores properties for accelerator creation.
#[derive(Clone)]
pub struct AcceleratorProps {
    /// Parameter set.
    pub params: ParamSet,

    /// Primitves.
    pub prims: Vec<ArcPrimitive>,
}
