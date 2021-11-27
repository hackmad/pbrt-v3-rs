//! Materials

#[macro_use]
extern crate lazy_static;

mod fourier;
mod glass;
mod matte;
mod mix;
mod plastic;

// Re-export
pub use fourier::*;
pub use glass::*;
pub use matte::*;
pub use mix::*;
pub use plastic::*;
