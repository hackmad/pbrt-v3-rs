//! Film

#![allow(dead_code)]
use super::geometry::*;
use super::pbrt::*;

/// Stores all of the sample values needed to specify a camera ray.
#[derive(Copy, Clone, Default)]
pub struct Film {
    /// The overall resolution of the image in pixels.
    pub full_resolution: Point2i,

    /// The diaogonal of the film's physical area in meters.
    pub diagonal: Float,
}

impl Film {
    /// Returns the actual extent of the film in the scene.
    pub fn get_physical_extent(&self) -> Bounds2f {
        let aspect = self.full_resolution.y as Float / self.full_resolution.x as Float;
        let x = (self.diagonal * self.diagonal / (1.0 + aspect * aspect)).sqrt();
        let y = aspect * x;
        Bounds2f::new(
            Point2f::new(-x / 2.0, -y / 2.0),
            Point2f::new(x / 2.0, y / 2.0),
        )
    }
}
