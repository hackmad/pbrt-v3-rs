//! Samplers

mod halton;
mod maxmin;
mod random;
mod sobol;
mod stratified;
mod zero_two_sequence;

// Re-export.
pub use halton::*;
pub use maxmin::*;
pub use random::*;
pub use sobol::*;
pub use stratified::*;
pub use zero_two_sequence::*;

use crate::core::geometry::Bounds2i;
use crate::core::paramset::ParamSet;

/// Stores properties for sampler creation.
#[derive(Clone)]
pub struct SamplerProps {
    /// Parameter set.
    pub params: ParamSet,

    /// Sample bounds.
    pub sample_bounds: Bounds2i,
}
