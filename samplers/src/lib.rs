//! Samplers

#[macro_use]
extern crate log;

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
