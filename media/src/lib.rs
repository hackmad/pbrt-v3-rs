//! Participating Media

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate log;

mod grid;
mod homogeneous;

// Re-export
pub use grid::*;
pub use homogeneous::*;
