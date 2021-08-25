//! Materials

#[macro_use]
extern crate lazy_static;

mod fourier;
mod matte;
mod mix;
mod plastic;

// Re-export
pub use fourier::*;
pub use matte::*;
pub use mix::*;
pub use plastic::*;
