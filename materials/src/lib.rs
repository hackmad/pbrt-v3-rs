//! Materials

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate log;

mod fourier;
mod glass;
mod kdsubsurface;
mod matte;
mod metal;
mod mirror;
mod mix;
mod plastic;
mod substrate;
mod subsurface;
mod translucent;
mod uber;

// Re-export
pub use fourier::*;
pub use glass::*;
pub use kdsubsurface::*;
pub use matte::*;
pub use metal::*;
pub use mirror::*;
pub use mix::*;
pub use plastic::*;
pub use substrate::*;
pub use subsurface::*;
pub use translucent::*;
pub use uber::*;
