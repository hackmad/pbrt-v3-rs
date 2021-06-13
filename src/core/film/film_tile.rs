//! Film tile

use crate::core::film::FILTER_TABLE_WIDTH;
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;

/// Stores contributions for the pixels in a region of the image.
pub struct FilmTile<'a> {
    /// Contributions of all pixels in the tile.
    pub pixels: Vec<FilmTilePixel>,

    /// Bounds of the pixels in the final image.
    pixel_bounds: Bounds2i,

    /// Filter radius.
    filter_radius: Vector2f,

    /// Reciprocal of `filter_radius`.
    inv_filter_radius: Vector2f,

    /// Filter table.
    filter_table: &'a [Float],

    /// Maximum sample luminence.
    max_sample_luminance: Float,
}

impl<'a> FilmTile<'a> {
    /// Create a new `FilmTile` instance.
    ///
    /// * `pixel_bounds`         - Bounds of the pixels in the final image.
    /// * `filter_radius`        - Filter radius.
    /// * `filter_table`         - Filter table.
    /// * `max_sample_luminance` - Optional maximum sample luminence to use use.
    ///                            Defaults to `INFINITY`.
    pub fn new(
        pixel_bounds: Bounds2i,
        filter_radius: Vector2f,
        filter_table: &'a [Float],
        max_sample_luminance: Option<Float>,
    ) -> Self {
        Self {
            pixel_bounds,
            filter_radius,
            inv_filter_radius: Vector2f::new(1.0 / filter_radius.x, 1.0 / filter_radius.y),
            filter_table,
            pixels: vec![FilmTilePixel::default(); max(0, pixel_bounds.area() as usize)],
            max_sample_luminance: match max_sample_luminance {
                Some(luminence) => luminence,
                None => INFINITY,
            },
        }
    }

    /// Add the radiance carried by a ray for a sample. This should be called by
    /// integrators.
    ///
    /// * `p_film`         - Point on film.
    /// * `l`              - Radiance value `L`.
    /// * `sample_weight`  - Weight for the sample's contribution.
    pub fn add_sample(&mut self, p_film: Point2f, l: Spectrum, sample_weight: Float) {
        let ly = l.y();
        let l = if ly > self.max_sample_luminance {
            l * self.max_sample_luminance / ly
        } else {
            l
        };

        // Compute sample's raster bounds.
        let p_film_discrete = p_film - Vector2f::new(0.5, 0.5);
        let mut p0 = Point2i::from((p_film_discrete - self.filter_radius).ceil());
        let mut p1 =
            Point2i::from((p_film_discrete + self.filter_radius).floor()) + Point2i::new(1, 1);
        p0 = p0.max(&self.pixel_bounds.p_min);
        p1 = p1.min(&self.pixel_bounds.p_max);

        // Loop over filter support and add sample to pixel arrays.
        let filter_table_size = FILTER_TABLE_WIDTH; // NOTE: not the entire size of the filter table.

        // Precompute `x` and `y` filter table offsets.
        let ifx: Vec<usize> = (p0.x..p1.x)
            .map(|x| {
                let fx = abs((x as Float - p_film_discrete.x)
                    * self.inv_filter_radius.x
                    * filter_table_size as Float);
                min(fx.floor(), filter_table_size as Float - 1.0) as usize
            })
            .collect();

        let ify: Vec<usize> = (p0.y..p1.y)
            .map(|y| {
                let fy = abs((y as Float - p_film_discrete.y)
                    * self.inv_filter_radius.y
                    * filter_table_size as Float);
                min(fy.floor(), filter_table_size as Float - 1.0) as usize
            })
            .collect();

        for y in p0.y..p1.y {
            for x in p0.x..p1.x {
                // Evaluate filter value at `(x, y)` pixel.
                let offset =
                    ify[(y - p0.y) as usize] * filter_table_size + ifx[(x - p0.x) as usize];
                let filter_weight = self.filter_table[offset];

                // Update pixel values with filtered sample contribution.
                let pixel_offset = self.get_pixel_offset(&Point2i::new(x, y));

                self.pixels[pixel_offset].contrib_sum += l * sample_weight * filter_weight;
                self.pixels[pixel_offset].filter_weight_sum += filter_weight;
            }
        }
    }

    /// Converts pixel coordinates with respect to the overall image and to
    /// coordinates in the film tile and returns the correspdoning pixel.
    ///
    /// * `p` - The pixel coordinates with respect to the overall image.
    pub fn get_pixel_offset(&self, p: &Point2i) -> usize {
        assert!(self.pixel_bounds.contains_exclusive(p));
        let width = self.pixel_bounds.p_max.x - self.pixel_bounds.p_min.x;
        let offset = (p.x - self.pixel_bounds.p_min.x) + (p.y - self.pixel_bounds.p_min.y) * width;
        offset as usize
    }

    /// Returns the bounds of the pixels in the final image.
    pub fn get_pixel_bounds(&self) -> Bounds2i {
        self.pixel_bounds
    }
}

/// Stores the contributions for pixel samples in a `FilmTile`.
#[derive(Default, Copy, Clone)]
pub struct FilmTilePixel {
    /// Sum of weighted contributions form the pixel samples.
    pub contrib_sum: Spectrum,

    /// Sum of filter weights.
    pub filter_weight_sum: Float,
}
