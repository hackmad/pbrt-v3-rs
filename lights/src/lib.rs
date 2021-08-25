//! Lights

#[macro_use]
extern crate log;

mod diffuse;
mod distant;
mod infinite;
mod point;

// Re-export.
pub use diffuse::*;
pub use distant::*;
pub use infinite::*;
pub use point::*;
