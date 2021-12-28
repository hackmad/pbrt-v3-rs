//! MIPMap

use super::MIPMap;
use crate::geometry::*;
use crate::image_io::write_image;
use crate::pbrt::*;
use crate::spectrum::RGBSpectrum;

/// Image writer interface for MIPMap levels.
pub trait MIPMapImageWriter {
    /// Write mipmap levels to images in the given base path.
    ///
    /// * `pyramids`  - The MIPMap pyramids.
    /// * `base_path` - Base path.
    fn write_images(&self, base_path: &str);
}

impl MIPMapImageWriter for MIPMap<RGBSpectrum> {
    /// Write mipmap levels to images in the given base path.
    ///
    /// * `pyramids`  - The MIPMap pyramids.
    /// * `base_path` - Base path.
    fn write_images(&self, base_path: &str) {
        for (i, level) in self.pyramid.iter().enumerate() {
            let nx = level.u_size();
            let ny = level.v_size();
            let n = 3 * nx * ny;

            let mut rgb = vec![0.0; n];
            let mut j = 0;
            for pixel in level.linear_vec() {
                rgb[j] = pixel[0];
                rgb[j + 1] = pixel[1];
                rgb[j + 2] = pixel[2];
                j += 3;
            }

            let path = format!("{}-{}.png", base_path, i);
            let bounds = Bounds2i::new(Point2i::new(0, 0), Point2i::new(nx as Int, ny as Int));
            write_image(&path, &rgb, &bounds).unwrap();
        }
    }
}

impl MIPMapImageWriter for MIPMap<Float> {
    /// Write mipmap levels to images in the given base path.
    ///
    /// * `pyramids`  - The MIPMap pyramids.
    /// * `base_path` - Base path.
    fn write_images(&self, base_path: &str) {
        for (i, level) in self.pyramid.iter().enumerate() {
            let nx = level.u_size();
            let ny = level.v_size();
            let n = 3 * nx * ny;

            let mut rgb = vec![0.0; n];
            let mut j = 0;
            for pixel in level.linear_vec() {
                rgb[j] = pixel;
                rgb[j + 1] = pixel;
                rgb[j + 2] = pixel;
                j += 3;
            }

            let path = format!("{}-{}.png", base_path, i);
            let bounds = Bounds2i::new(Point2i::new(0, 0), Point2i::new(nx as Int, ny as Int));
            write_image(&path, &rgb, &bounds).unwrap();
        }
    }
}
