//! Film

use crate::app::OPTIONS;
use crate::app::WINDOW_HEIGHT;
use crate::app::WINDOW_PIXELS;
use crate::app::WINDOW_WIDTH;
use crate::filter::*;
use crate::geometry::*;
use crate::image_io::*;
use crate::paramset::*;
use crate::pbrt::*;
use crate::spectrum::*;
use crate::{stat_inc, stat_memory_counter, stat_register_fns, stats::*};
use std::sync::RwLockWriteGuard;
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

stat_memory_counter!("Memory/Film pixels", FILM_PIXEL_MEMORY, film_stats_pixels);

stat_register_fns!(film_stats_pixels);

/// Pixel data.
#[derive(Copy, Clone, Default)]
#[repr(C)]
pub struct Pixel {
    /// Stores the running weighted sums of spectral pixel contributions using XYZ colors.
    pub xyz: [Float; 3],

    /// Holds the sum of filter weight values for the sample contributions to the pixel.
    pub filter_weight_sum: Float,

    /// Holds an unweighted sum of sample splats.
    pub splat_xyz: [Float; 3],

    /// Used to pad this struct to the 32-bit/64-bit. This will work for both `Float` => `f32` and `Float` => `f64`.
    pad: Float,
}

/// Models the sensing device in a simulated camera. It stores all of the sample values needed to specify a camera ray.
pub struct Film {
    /// The overall image resolution in pixels.
    pub full_resolution: Point2i,

    /// The diaogonal of the film's physical area in meters.
    pub diagonal: Float,

    /// Filter function to use for image reconstruction from samples.
    pub filter: ArcFilter,

    /// Filename of output image.
    pub filename: String,

    /// Crop window of the subset of the image to render.
    pub cropped_pixel_bounds: Bounds2i,

    /// The filter table.
    filter_table: Arc<[Float; FILTER_TABLE_SIZE]>,

    /// Scale factor for pixel values.
    scale: Float,

    /// Maximum sample luminence.
    max_sample_luminance: Float,

    /// Stores the image pixels.
    pixels: RwLock<Vec<Pixel>>,
}

impl Film {
    /// Create a new `Film` instance.
    ///
    /// * `resolution`           - The overall image resolution in pixels.
    /// * `crop_window`          - Crop window of the subset of the image to render.
    /// * `filter`               - Filter function to use for image reconstruction from samples.
    /// * `diagonal`             - The diaogonal of the film's physical area in millimeters.
    /// * `filename`             - Filename of output image.
    /// * `scale`                - Optional scale factor for pixel values. If None specified, sets to 1.0.
    /// * `max_sample_luminance` - Optional maximum sample luminence to use use. Defaults to `INFINITY`.
    pub fn new(
        resolution: &Point2i,
        crop_window: &Bounds2f,
        filter: ArcFilter,
        diagonal: Float,
        filename: &str,
        scale: Option<Float>,
        max_sample_luminance: Option<Float>,
    ) -> Self {
        register_stats();

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
        let filter_data = filter.get_data();
        let mut filter_table = [0.0; FILTER_TABLE_SIZE];
        let mut offset = 0;
        for y in 0..FILTER_TABLE_WIDTH {
            for x in 0..FILTER_TABLE_WIDTH {
                let p = Point2f::new(
                    (x as Float + 0.5) * filter_data.radius.x * INV_FILTER_TABLE_WIDTH,
                    (y as Float + 0.5) * filter_data.radius.y * INV_FILTER_TABLE_WIDTH,
                );
                filter_table[offset] = filter.evaluate(&p);
                offset += 1;
            }
        }

        // Allocate film image storage.
        let n = cropped_pixel_bounds.area() as usize;
        let pixels = RwLock::new(vec![Pixel::default(); n]);
        stat_inc!(FILM_PIXEL_MEMORY, (n * std::mem::size_of::<Pixel>()) as u64);

        Self {
            full_resolution: *resolution,
            diagonal: diagonal * 0.001, // Convert to meters.
            filter,
            filter_table: Arc::new(filter_table),
            filename: String::from(filename),
            cropped_pixel_bounds,
            scale: scale.unwrap_or(1.0),
            max_sample_luminance: match max_sample_luminance {
                Some(luminence) => luminence,
                None => INFINITY,
            },
            pixels,
        }
    }

    /// Returns the sample bounds accounting for the half-pixel offsets when converting from discrete to continuous
    /// pixel coordinates.
    pub fn get_sample_bounds(&self) -> Bounds2i {
        let filter_data = self.filter.get_data();
        let half_pixel = Vector2f::new(0.5, 0.5);

        let p0 = (Point2f::from(self.cropped_pixel_bounds.p_min) + half_pixel - filter_data.radius).floor();
        let p1 = (Point2f::from(self.cropped_pixel_bounds.p_max) - half_pixel + filter_data.radius).ceil();
        let float_bounds = Bounds2f::new(p0, p1);

        Bounds2i::from(float_bounds)
    }

    /// Returns the actual extent of the film in the scene.
    pub fn get_physical_extent(&self) -> Bounds2f {
        let aspect = self.full_resolution.y as Float / self.full_resolution.x as Float;
        let x = (self.diagonal * self.diagonal / (1.0 + aspect * aspect)).sqrt();
        let y = aspect * x;
        Bounds2f::new(Point2f::new(-x / 2.0, -y / 2.0), Point2f::new(x / 2.0, y / 2.0))
    }

    /// Gets the pixel given its coordinates in the overall image.
    ///
    /// * `p` - The pixel coordinates with respect to the overall image.
    pub fn get_pixel_offset(&self, p: &Point2i) -> usize {
        assert!(self.cropped_pixel_bounds.contains_exclusive(p));
        let width = self.cropped_pixel_bounds.p_max.x - self.cropped_pixel_bounds.p_min.x;
        let offset = (p.x - self.cropped_pixel_bounds.p_min.x) + (p.y - self.cropped_pixel_bounds.p_min.y) * width;
        offset as usize
    }

    /// Returns a `FilmTile` that stores the contributions for pixels in the specified region of the image.
    ///
    /// * `sample_bounds` - Tile region in the overall image.
    pub fn get_film_tile(&self, sample_bounds: Bounds2i) -> FilmTile {
        let filter_data = self.filter.get_data();
        let half_pixel = Vector2f::new(0.5, 0.5);

        // Bound image pixels that samples in `sample_bounds` contribute to.
        let float_bounds = Bounds2f::from(sample_bounds);
        let p0 = Point2i::from((float_bounds.p_min - half_pixel - filter_data.radius).ceil());
        let p1 = Point2i::from((float_bounds.p_max - half_pixel + filter_data.radius).floor()) + Point2i::new(1, 1);
        let tile_pixel_bounds = Bounds2i::new(p0, p1).intersect(&self.cropped_pixel_bounds);

        FilmTile::new(
            tile_pixel_bounds,
            filter_data.radius,
            Arc::clone(&self.filter_table),
            Some(self.max_sample_luminance),
        )
    }

    /// Clear the splats for all pixels in the image.
    pub fn clear(&self) {
        let mut pixels = self.pixels.write().unwrap();
        for pixel in self.cropped_pixel_bounds {
            let pixel_offset = self.get_pixel_offset(&pixel);
            pixels[pixel_offset].splat_xyz = [0.0; 3];
            pixels[pixel_offset].filter_weight_sum = 0.0;
        }
    }

    /// Merge the `FilmTile`'s pixel contribution into the image. Also merge it with the render preview in window.
    ///
    /// This is same as `merge_film_tile_old()` to reduce all but one call to `tile.get_pixel_offset()` and
    /// `self.get_pixel_offset()`.
    ///
    /// * `tile`        - The `FilmTile` to merge.
    /// * `splat_scale` - Scale factor for `add_splat()` (default = 1.0).
    pub fn merge_film_tile(&self, tile: &FilmTile, splat_scale: Float) {
        let mut pixels = self.pixels.write().unwrap();

        let cropped_pixel_width = (self.cropped_pixel_bounds.p_max.x - self.cropped_pixel_bounds.p_min.x) as usize;

        let tile_pixel_bounds = tile.get_pixel_bounds();
        let tile_width = (tile_pixel_bounds.p_max.x - tile_pixel_bounds.p_min.x) as usize;
        let tile_size = tile_pixel_bounds.area() as usize;

        let mut tile_pixel = tile.get_pixel_offset(&tile_pixel_bounds.p_min);
        let mut merge_pixel = self.get_pixel_offset(&tile_pixel_bounds.p_min);
        let mut x = tile_pixel_bounds.p_min.x;

        for _ in 0..tile_size {
            let xyz = tile.pixels[tile_pixel].contrib_sum.to_xyz();
            pixels[merge_pixel].xyz[0] += xyz[0];
            pixels[merge_pixel].xyz[1] += xyz[1];
            pixels[merge_pixel].xyz[2] += xyz[2];

            pixels[merge_pixel].filter_weight_sum += tile.pixels[tile_pixel].filter_weight_sum;

            x += 1;
            tile_pixel += 1;
            merge_pixel += 1;

            if x == tile_pixel_bounds.p_max.x {
                x = tile_pixel_bounds.p_min.x;
                merge_pixel += cropped_pixel_width - tile_width;
            }
        }

        self.merge_film_tile_for_preview(
            tile_pixel_bounds,
            tile_size,
            tile_width,
            cropped_pixel_width,
            splat_scale,
            &mut pixels,
        );
    }

    /// Merge the `FilmTile`'s pixel contribution into the image.
    ///
    /// * `tile` - The `FilmTile` to merge.
    #[allow(unused)]
    fn merge_film_tile_old(&self, tile: &FilmTile) {
        let mut pixels = self.pixels.write().unwrap();
        for pixel in tile.get_pixel_bounds() {
            let tile_pixel = tile.get_pixel_offset(&pixel);
            let merge_pixel = self.get_pixel_offset(&pixel);
            let xyz = tile.pixels[tile_pixel].contrib_sum.to_xyz();
            for (i, colour) in xyz.iter().enumerate() {
                pixels[merge_pixel].xyz[i] += *colour;
            }
            pixels[merge_pixel].filter_weight_sum += tile.pixels[tile_pixel].filter_weight_sum;
        }
    }

    /// Merge the `FilmTile`'s pixel contribution into the window image for preview.
    ///
    /// * `tile`                - The `FilmTile` to merge.
    /// * `tile_size`           - Area of tile in pixels.
    /// * `tile_width`          - Horizontal size of tile.
    /// * `cropped_pixel_width` - Cropped width of the subset of the image to render.
    /// * `splat_scale`         - Scale factor for `add_splat()` (default = 1.0).
    /// * `pixels`              - The RwLockWriteGuard for `self.pixels`.
    fn merge_film_tile_for_preview(
        &self,
        tile_pixel_bounds: Bounds2i,
        tile_size: usize,
        tile_width: usize,
        cropped_pixel_width: usize,
        splat_scale: Float,
        pixels: &mut RwLockWriteGuard<Vec<Pixel>>,
    ) {
        let mut merge_pixel = self.get_pixel_offset(&tile_pixel_bounds.p_min);
        let mut x = tile_pixel_bounds.p_min.x;
        WINDOW_PIXELS
            .get()
            .map(|wp| wp.write().ok())
            .flatten()
            .map(|mut window_pixels| {
                let wp = window_pixels.frame_mut();
                for _ in 0..tile_size {
                    let pixel_rgb = self.get_pixel_rgb(
                        &pixels[merge_pixel].xyz,
                        &pixels[merge_pixel].splat_xyz,
                        pixels[merge_pixel].filter_weight_sum,
                        splat_scale,
                    );

                    // TODO Account for window dimensions not matching rendered image dimensions by scaling the
                    // tile appropriately.
                    let window_pixels_x = merge_pixel % WINDOW_WIDTH as usize;
                    let window_pixels_y = merge_pixel / WINDOW_HEIGHT as usize;
                    let offset = (window_pixels_y * WINDOW_WIDTH as usize + window_pixels_x) * 4; // RGBA

                    let rgb = apply_gamma(&pixel_rgb);
                    wp[offset + 0] = rgb[0];
                    wp[offset + 1] = rgb[1];
                    wp[offset + 2] = rgb[2];
                    wp[offset + 3] = 255; // Alpha

                    x += 1;
                    merge_pixel += 1;

                    if x == tile_pixel_bounds.p_max.x {
                        x = tile_pixel_bounds.p_min.x;
                        merge_pixel += cropped_pixel_width - tile_width;
                    }
                }
            });
    }

    /// Sets all pixel values in the cropped area with the given spectrum values.
    ///
    /// * `img` - The spectrum values for the cropped area.
    pub fn set_image(&self, img: &[Spectrum]) {
        let mut pixels = self.pixels.write().unwrap();
        let n_pixels = self.cropped_pixel_bounds.area() as usize;
        for i in 0..n_pixels {
            pixels[i].xyz = img[i].to_xyz();
            pixels[i].filter_weight_sum = 1.0;
            pixels[i].splat_xyz = [0.0; 3];
        }
    }

    /// Add `splat` contributions to a pixel.
    ///
    /// * `p` - The pixel coordinates with respect to the overall image.
    /// * `v` - `Splat` contribution to add to the pixel.
    pub fn add_splat(&self, p: &Point2f, v: &Spectrum) {
        if v.has_nans() {
            warn!("Ignoring splatted spectrum with NaN values at ({}, {})", p.x, p.y);
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
            let mut pixels = self.pixels.write().unwrap();
            for (i, colour) in xyz.iter().enumerate() {
                pixels[pixel_offset].splat_xyz[i] += colour;
            }
        }
    }

    /// Write the image to an output file.
    ///
    /// * `splat_scale` - Scale factor for `add_splat()` (default = 1.0).
    pub fn write_image(&self, splat_scale: Float) {
        info!("Converting image to RGB and computing final weighted pixel values");

        let n = 3 * self.cropped_pixel_bounds.area() as usize;
        let mut rgb = vec![0.0; n];

        let pixels = self.pixels.read().unwrap();
        for p in self.cropped_pixel_bounds {
            // Convert pixel XYZ color to RGB.
            let pixel_offset = self.get_pixel_offset(&p);

            let pixel_rgb = self.get_pixel_rgb(
                &pixels[pixel_offset].xyz,
                &pixels[pixel_offset].splat_xyz,
                pixels[pixel_offset].filter_weight_sum,
                splat_scale,
            );

            let rgb_offset = 3 * pixel_offset;
            rgb[rgb_offset + 0] = pixel_rgb[0];
            rgb[rgb_offset + 1] = pixel_rgb[1];
            rgb[rgb_offset + 2] = pixel_rgb[2];
        }

        // Write RGB image.
        if let Err(err) = write_image(&self.filename, &rgb, &self.cropped_pixel_bounds) {
            panic!("Error writing output image {}. {:}.", self.filename, err);
        }
    }

    /// Compute the RGB colour value for a pixel from XYZ colour space.
    ///
    /// * `pixel_xyz`         - The XYZ colour value of pixel.
    /// * `splat_xyz`         - The unweighted sum of sample splats.
    /// * `filter_weight_sum` - The sum of filter weight values for the sample contributions to the pixel.
    /// * `splat_scale`       - Scale factor for `add_splat()` (default = 1.0).
    fn get_pixel_rgb(
        &self,
        pixel_xyz: &[Float; 3],
        splat_xyz: &[Float; 3],
        filter_weight_sum: Float,
        splat_scale: Float,
    ) -> [Float; 3] {
        let mut rgb = xyz_to_rgb(pixel_xyz);
        let splat_rgb = xyz_to_rgb(splat_xyz);

        for v in rgb.iter_mut() {
            if filter_weight_sum != 0.0 {
                // Normalize pixel with weight sum.
                let inv_wt = 1.0 / filter_weight_sum;
                *v = max(0.0, *v * inv_wt);
            }

            // Add splat value at pixel.
            *v += splat_scale * splat_rgb[0];

            // Scale pixel value by `scale`.
            *v *= self.scale;
        }

        rgb
    }
}

impl From<(&ParamSet, ArcFilter)> for Film {
    /// Create a `BVHAccel` from given parameter set and filter.
    ///
    /// * `p` - Tuple containing the parameter set and filter.
    fn from(p: (&ParamSet, ArcFilter)) -> Self {
        let (params, filter) = p;

        let filename = if let Some(image_file) = OPTIONS.image_file.as_ref() {
            let params_filename = params.find_one_string("filename", "".to_owned());
            if !params_filename.is_empty() {
                warn!(
                    "Output filename supplied on command line, '{}' is overriding filename provided in scene description file, '{}'.",
                    image_file, params_filename
                );
            }
            image_file.clone()
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
        let mut crop = Bounds2f::EMPTY;
        if cwi == 4 {
            crop.p_min.x = clamp(min(cr[0], cr[1]), 0.0, 1.0);
            crop.p_max.x = clamp(max(cr[0], cr[1]), 0.0, 1.0);
            crop.p_min.y = clamp(min(cr[2], cr[3]), 0.0, 1.0);
            crop.p_max.y = clamp(max(cr[2], cr[3]), 0.0, 1.0);
        } else if cwi > 0 {
            panic!("{} values supplied for 'cropwindow'. Expected 4.", cwi);
        } else if OPTIONS.crop_window.len() == 4 {
            crop = Bounds2f::new(
                Point2f::new(
                    clamp(OPTIONS.crop_window[0], 0.0, 1.0),
                    clamp(OPTIONS.crop_window[2], 0.0, 1.0),
                ),
                Point2f::new(
                    clamp(OPTIONS.crop_window[1], 0.0, 1.0),
                    clamp(OPTIONS.crop_window[3], 0.0, 1.0),
                ),
            );
        } else {
            crop = Bounds2f::new(Point2f::new(0.0, 0.0), Point2f::new(1.0, 1.0));
        }

        let scale = params.find_one_float("scale", 1.0);
        let diagonal = params.find_one_float("diagonal", 35.0);
        let max_sample_luminance = params.find_one_float("maxsampleluminance", INFINITY);
        Self::new(
            &Point2i::new(xres, yres),
            &crop,
            filter,
            diagonal,
            &filename,
            Some(scale),
            Some(max_sample_luminance),
        )
    }
}
