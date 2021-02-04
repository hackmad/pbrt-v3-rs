//! Film

#![allow(dead_code)]
use crate::core::app::OPTIONS;
use crate::core::filter::*;
use crate::core::geometry::*;
use crate::core::image_io::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;
use std::sync::{Arc, RwLock};

mod film_tile;

// Re-export.
pub use film_tile::*;

/// Filter table width.
pub const FILTER_TABLE_WIDTH: usize = 16;

/// Filter table size.
pub const FILTER_TABLE_SIZE: usize = FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH;

/// Reciprocal of `FILTER_TABLE_WIDTH`.
pub const INV_FILTER_TABLE_WIDTH: Float = 1.0 / (FILTER_TABLE_WIDTH as Float);

/// Pixel data.
#[derive(Copy, Clone, Default)]
#[repr(C)]
pub struct Pixel {
    /// Stores the running weighted sums of spectral pixel contributions using
    /// XYZ colors.
    pub xyz: [Float; 3],

    /// Holds the sum of filter weight values for the sample contributions to
    /// the pixel.
    pub filter_weight_sum: Float,

    /// Holds an unweighted sum of sample splats.
    pub splat_xyz: [Float; 3],

    /// Used to pad this struct to the 32-bit/64-bit. This will work for both
    /// `Float` => `f32` and `Float` => `f64`.
    pad: Float,
}

/// Models the sensing device in a simulated camera. It stores all of the sample
/// values needed to specify a camera ray.
#[derive(Clone)]
pub struct Film {
    /// The overall image resolution in pixels.
    pub full_resolution: Point2i,

    /// The diaogonal of the film's physical area in meters.
    pub diagonal: Float,

    /// Filter function to use for image reconstruction from samples.
    pub filter: ArcFilter,

    /// Filename of output image.
    pub filename: String,

    /// Crop windo of the subset of the image to render.
    pub cropped_pixel_bounds: Bounds2i,

    /// The filter table.
    filter_table: [Float; FILTER_TABLE_SIZE],

    /// Scale factor for pixel values.
    scale: Float,

    /// Maximum sample luminence.
    max_sample_luminance: Float,

    /// Stores the image pixels.
    pixels: Arc<RwLock<Vec<Pixel>>>,
}

impl Film {
    /// Create a new `Film` instance.
    ///
    /// * `resolution`           - The overall image resolution in pixels.
    /// * `crop_window`          - Crop window of the subset of the image to render.
    /// * `filter`               - Filter function to use for image reconstruction
    ///                            from samples.
    /// * `diagonal`             - The diaogonal of the film's physical area in
    //                             millimeters.
    /// * `filename`             - Filename of output image.
    /// * `scale`                - Optional scale factor for pixel values. If
    ///                            None specified, sets to 1.0.
    /// * `max_sample_luminance` - Optional maximum sample luminence to use use.
    ///                            Defaults to `INFINITY`.
    pub fn new(
        resolution: &Point2i,
        crop_window: &Bounds2f,
        filter: ArcFilter,
        diagonal: Float,
        filename: &str,
        scale: Option<Float>,
        max_sample_luminance: Option<Float>,
    ) -> Self {
        // Compute the film image bounds.
        let cropped_pixel_bounds = Bounds2i::new(
            Point2i::new(
                (resolution.x as Float * crop_window.p_min.x).ceil() as Int,
                (resolution.y as Float * crop_window.p_min.y).ceil() as Int,
            ),
            Point2i::new(
                (resolution.x as Float * crop_window.p_max.x).ceil() as Int,
                (resolution.y as Float * crop_window.p_max.y).ceil() as Int,
            ),
        );

        // Precompute filter weight table.
        let filter = filter.clone();
        let filter_data = filter.get_data();
        let mut filter_table = [0.0; FILTER_TABLE_SIZE];
        let mut offset = 0;
        for y in 0..FILTER_TABLE_WIDTH {
            for x in 0..FILTER_TABLE_WIDTH {
                let p = Point2f::new(
                    (x as Float + 0.5) * filter_data.radius.x / INV_FILTER_TABLE_WIDTH,
                    (y as Float + 0.5) * filter_data.radius.y / INV_FILTER_TABLE_WIDTH,
                );
                filter_table[offset] = filter.evaluate(&p);
                offset += 1;
            }
        }

        // Allocate film image storage.
        let n = cropped_pixel_bounds.area() as usize;
        let pixels = Arc::new(RwLock::new(vec![Pixel::default(); n]));

        Self {
            full_resolution: *resolution,
            diagonal: diagonal * 0.001, // Convert to meters.
            filter,
            filter_table,
            filename: String::from(filename),
            cropped_pixel_bounds,
            scale: match scale {
                Some(s) => s,
                None => 1.0,
            },
            max_sample_luminance: match max_sample_luminance {
                Some(luminence) => luminence,
                None => INFINITY,
            },
            pixels,
        }
    }

    /// Returns the sample bounds accounting for the half-pixel offsets when
    /// converting from discrete to continuous pixel coordinates.
    pub fn get_sample_bounds(&self) -> Bounds2i {
        let filter_data = self.filter.get_data();
        let half_pixel = Vector2f::new(0.5, 0.5);

        let p0 = (Point2f::from(self.cropped_pixel_bounds.p_min) + half_pixel - filter_data.radius)
            .floor();
        let p1 = (Point2f::from(self.cropped_pixel_bounds.p_max) - half_pixel + filter_data.radius)
            .ceil();
        let float_bounds = Bounds2f::new(p0, p1);

        Bounds2i::from(float_bounds)
    }

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

    /// Gets the pixel given its coordinates in the overall image.
    ///
    /// * `p` - The pixel coordinates with respect to the overall image.
    pub fn get_pixel_offset(&self, p: &Point2i) -> usize {
        assert!(self.cropped_pixel_bounds.contains_exclusive(p));
        let width = self.cropped_pixel_bounds.p_max.x - self.cropped_pixel_bounds.p_min.x;
        let offset = (p.x - self.cropped_pixel_bounds.p_min.x)
            + (p.y - self.cropped_pixel_bounds.p_min.y) * width;
        offset as usize
    }

    /// Returns an `Arc<FilmTile>` that stores the contributions for pixels in
    /// the specified region of the image.
    ///
    /// * `sample_bounds` - Tile region in the overall image.
    pub fn get_film_tile(&self, sample_bounds: Bounds2i) -> Arc<FilmTile> {
        let filter_data = self.filter.get_data();
        let half_pixel = Vector2f::new(0.5, 0.5);

        // Bound image pixels that samples in `sample_bounds` contribute to.
        let float_bounds = Bounds2f::from(sample_bounds);
        let p0 = Point2i::from((float_bounds.p_min - half_pixel - filter_data.radius).ceil());
        let p1 = Point2i::from((float_bounds.p_max - half_pixel + filter_data.radius).floor())
            + Point2i::new(1, 1);
        let tile_pixel_bounds = Bounds2i::new(p0, p1).intersect(&self.cropped_pixel_bounds);

        Arc::new(FilmTile::new(
            tile_pixel_bounds,
            filter_data.radius,
            &self.filter_table,
            Some(self.max_sample_luminance),
        ))
    }

    /// Clear the splats for all pixels in the image.
    pub fn clear(&mut self) {
        let mut pixels = self.pixels.write().unwrap();
        for pixel in self.cropped_pixel_bounds {
            let pixel_offset = self.get_pixel_offset(&pixel);
            (*pixels)[pixel_offset].splat_xyz = [0.0; 3];
        }
    }

    /// Merge the `FilmTile`'s pixel contribution into the image.
    ///
    /// * `tile` - The `FilmTile` to merge.
    pub fn merge_film_tile(&mut self, tile: Arc<FilmTile>) {
        let mut pixels = self.pixels.write().unwrap();
        for pixel in tile.get_pixel_bounds() {
            let tile_pixel = tile.get_pixel_offset(&pixel);
            let merge_pixel = self.get_pixel_offset(&pixel);
            let xyz = tile.pixels[tile_pixel].contrib_sum.to_xyz();
            for i in 0..3 {
                (*pixels)[merge_pixel].xyz[i] += xyz[i];
            }
            (*pixels)[merge_pixel].filter_weight_sum += tile.pixels[tile_pixel].filter_weight_sum;
        }
    }

    /// Sets all pixel values in the cropped area with the given spectrum values.
    ///
    /// * `img` - The spectrum values for the cropped area.
    pub fn set_image(&mut self, img: &[Spectrum]) {
        let mut pixels = self.pixels.write().unwrap();
        let n_pixels = self.cropped_pixel_bounds.area();
        for i in (0..n_pixels).map(|i| i as usize) {
            (*pixels)[i].xyz = img[i].to_xyz();
            (*pixels)[i].filter_weight_sum = 1.0;
            (*pixels)[i].splat_xyz = [0.0; 3];
        }
    }

    /// Add `splat` contributions to a pixel.
    ///
    /// * `p` - The pixel coordinates with respect to the overall image.
    /// * `v` - `Splat` contribution to add to the pixel.
    pub fn add_splat(&mut self, p: &Point2f, v: &Spectrum) {
        if v.has_nans() {
            warn!(
                "Ignoring splatted spectrum with NaN values at ({}, {})",
                p.x, p.y
            );
            return;
        }

        let vy = v.y();
        if vy < 0.0 {
            warn!(
                "Ignoring splatted spectrum with negative luminance {} at ({}, {})",
                vy, p.x, p.y
            );
        } else if vy.is_infinite() {
            warn!(
                "Ignoring splatted spectrum with infinite luminance at ({}, {})",
                p.x, p.y
            );
        } else {
            let mut pixels = self.pixels.write().unwrap();

            let pi = Point2i::from(p.floor());
            if !self.cropped_pixel_bounds.contains_exclusive(&pi) {
                return;
            }

            let v = if vy > self.max_sample_luminance {
                *v * self.max_sample_luminance / vy
            } else {
                *v
            };

            let xyz = v.to_xyz();
            let pixel_offset = self.get_pixel_offset(&pi);
            for i in 0..3 {
                (*pixels)[pixel_offset].splat_xyz[i] += xyz[i];
            }
        }
    }

    /// Write the image to an output file.
    ///
    /// * `splat_scale` - Scale factor provided to `add_splat()`.
    pub fn write_image(&self, splat_scale: Float) {
        info!("Converting image to RGB and computing final weighted pixel values");

        let mut pixels = self.pixels.write().unwrap();

        let n = 3 * self.cropped_pixel_bounds.area() as usize;
        let mut rgb = vec![0.0; n];

        let mut offset = 0;
        for p in self.cropped_pixel_bounds {
            // Convert pixel XYZ color to RGB.
            let pixel_offset = self.get_pixel_offset(&p);
            (*pixels)[pixel_offset].xyz =
                xyz_to_rgb(&[rgb[3 * offset], rgb[3 * offset + 1], rgb[3 * offset + 2]]);

            // Normalize pixel with weight sum.
            let filter_weight_sum = (*pixels)[pixel_offset].filter_weight_sum;
            if filter_weight_sum != 0.0 {
                let inv_wt = 1.0 / filter_weight_sum;
                rgb[3 * offset] = max(0.0, rgb[3 * offset] * inv_wt);
                rgb[3 * offset + 1] = max(0.0, rgb[3 * offset + 1] * inv_wt);
                rgb[3 * offset + 2] = max(0.0, rgb[3 * offset + 2] * inv_wt);
            }

            // Add splat value at pixel.
            let splat_rgb = xyz_to_rgb(&(*pixels)[pixel_offset].splat_xyz);
            rgb[3 * offset] += splat_scale * splat_rgb[0];
            rgb[3 * offset + 1] += splat_scale * splat_rgb[1];
            rgb[3 * offset + 2] += splat_scale * splat_rgb[2];

            // Scale pixel value by `scale`.
            rgb[3 * offset] *= self.scale;
            rgb[3 * offset + 1] *= self.scale;
            rgb[3 * offset + 2] *= self.scale;

            offset += 1;
        }

        // Write RGB image
        if let Err(err) = write_image(&self.filename, &rgb, &self.cropped_pixel_bounds) {
            panic!("Error writing output image {}. {:}.", self.filename, err);
        }
    }
}

impl From<(&mut ParamSet, ArcFilter)> for Film {
    /// Create a `BVHAccel` from given parameter set and filter.
    ///
    /// * `p` - Tuple containing the parameter set and filter.
    fn from(p: (&mut ParamSet, ArcFilter)) -> Self {
        let (params, filter) = p;

        let image_file = &OPTIONS.image_file[..];
        let filename = if image_file.len() > 0 {
            let params_filename = params.find_one_string("filename", String::from(""));
            if params_filename.len() > 0 {
                warn!(
                    "Output filename supplied on command line, '{}' is overriding 
                    filename provided in scene description file, '{}'.",
                    OPTIONS.image_file, params_filename
                );
                params_filename
            } else {
                String::from(image_file)
            }
        } else {
            params.find_one_string("filename", String::from("pbrt.exr"))
        };

        let mut xres = params.find_one_int("xresolution", 1280);
        let mut yres = params.find_one_int("yresolution", 720);
        if OPTIONS.quick_render {
            xres = max(1, xres / 4);
            yres = max(1, yres / 4);
        }

        let cr = params.find_float("cropwindow");
        let cwi = cr.len();
        let mut crop = Bounds2f::default();
        if cwi == 4 {
            crop.p_min.x = clamp(min(cr[0], cr[1]), 0.0, 1.0);
            crop.p_max.x = clamp(max(cr[0], cr[1]), 0.0, 1.0);
            crop.p_min.y = clamp(min(cr[2], cr[3]), 0.0, 1.0);
            crop.p_max.y = clamp(max(cr[2], cr[3]), 0.0, 1.0);
        } else if cwi > 0 {
            panic!("{} values supplied for 'cropwindow'. Expected 4.", cwi);
        } else {
            crop = Bounds2f::new(
                Point2f::new(
                    clamp(OPTIONS.crop_window[0][0], 0.0, 1.0),
                    clamp(OPTIONS.crop_window[1][0], 0.0, 1.0),
                ),
                Point2f::new(
                    clamp(OPTIONS.crop_window[0][1], 0.0, 1.0),
                    clamp(OPTIONS.crop_window[1][1], 0.0, 1.0),
                ),
            );
        }

        let scale = params.find_one_float("scale", 1.0);
        let diagonal = params.find_one_float("diagonal", 35.0);
        let max_sample_luminance = params.find_one_float("maxsampleluminance", INFINITY);
        Self::new(
            &Point2i::new(xres, yres),
            &crop,
            filter.clone(),
            diagonal,
            &filename,
            Some(scale),
            Some(max_sample_luminance),
        )
    }
}
