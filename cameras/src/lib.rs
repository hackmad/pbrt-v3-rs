//! Camera

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate log;

mod environment_camera;
mod orthographic_camera;
mod perspective_camera;
mod realistic_camera;

// Re-export
pub use environment_camera::*;
pub use orthographic_camera::*;
pub use perspective_camera::*;
pub use realistic_camera::*;
