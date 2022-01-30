//! Lights

#[macro_use]
extern crate log;

mod diffuse;
mod distant;
mod goniometric;
mod infinite;
mod point;
mod projection;
mod spot;

// Re-export.
pub use diffuse::*;
pub use distant::*;
pub use goniometric::*;
pub use infinite::*;
pub use point::*;
pub use projection::*;
pub use spot::*;
