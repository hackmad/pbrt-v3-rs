//! Materials

#[macro_use]
extern crate lazy_static;

mod fourier;
mod glass;
mod matte;
mod metal;
mod mirror;
mod mix;
mod plastic;
mod translucent;
mod uber;

// Re-export
pub use fourier::*;
pub use glass::*;
pub use matte::*;
pub use metal::*;
pub use mirror::*;
pub use mix::*;
pub use plastic::*;
pub use translucent::*;
pub use uber::*;
