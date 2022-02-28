//! Integrators

#[macro_use]
extern crate log;

mod direct_lighting;
mod whitted;

// Re-export.
pub use direct_lighting::*;
pub use whitted::*;
