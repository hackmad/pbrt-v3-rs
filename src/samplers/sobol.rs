//! Sobol Sampler.

use crate::core::app::OPTIONS;
use crate::core::geometry::*;
use crate::core::low_discrepency::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;
use crate::core::rng::*;
use crate::core::sampler::*;
use crate::core::sobol_matrices::*;
use std::sync::Arc;

/// Implements the Sobol Sampler.
pub struct SobolSampler {
    /// The sampler data.
    data: SamplerData,

    /// The global sampler data.
    gdata: GlobalSamplerData,

    /// sample bounds.
    sample_bounds: Bounds2i,

    /// Stores the next power of two based on maximum extend in x or y direction
    /// of `sample_bounds`.
    resolution: i32,

    /// Log base 2 of `resolution`.
    log_2_resolution: i32,
}

impl SobolSampler {
    /// Create a new `SobolSampler`.
    ///
    /// * `samples_per_pixel` - Number of samples per pixel.
    /// * `sample_bounds`     - Sample bounds.
    fn new(samples_per_pixel: usize, sample_bounds: Bounds2i) -> Self {
        let resolution = max(sample_bounds.diagonal().x, sample_bounds.diagonal().y);

        Self {
            data: SamplerData::new(samples_per_pixel),
            gdata: GlobalSamplerData::new(),
            sample_bounds,
            resolution,
            log_2_resolution: resolution.log2(),
        }
    }
}

impl SobolSampler {
    /// Performs the inverse mapping from the current pixel and given sample
    /// index to a global index into the overall set of sample vectors.
    ///
    /// * `sample_num` - The sample number.
    fn get_index_for_sample(&mut self, sample_num: usize) -> u64 {
        sobol_interval_to_index(
            self.log_2_resolution as u32,
            sample_num as u64,
            &Point2i::from(self.data.current_pixel - self.sample_bounds.p_min),
        )
    }

    /// Returns the sample value for the given dimension of the index^th sample
    /// vector in the sequence.
    ///
    /// * `index` - Index of the sample.
    /// * `dim`   - Dimension.
    fn sample_dimension(&mut self, index: u64, dim: u16) -> Float {
        assert!(
            (dim as usize) <= NUM_SOBOL_DIMENSIONS,
            "SobolSampler can only sample up to {} dimensions.",
            NUM_SOBOL_DIMENSIONS
        );

        let mut s = sobol_sample(index, dim, 0);
        if dim == 0 || dim == 1 {
            s = s * (self.resolution as Float) + self.sample_bounds.p_min[dim as usize] as Float;
            s = clamp(
                s - self.data.current_pixel[dim as usize] as Float,
                0.0,
                ONE_MINUS_EPSILON,
            );
        }
        s
    }
}

impl Sampler for SobolSampler {
    /// Returns the underlying `SamplerData`.
    fn get_data(&mut self) -> &mut SamplerData {
        &mut self.data
    }

    /// Generates a new instance of an initial `Sampler` for use by a rendering
    /// thread.
    ///
    /// * `seed` - The seed for the random number generator (ignored).
    fn clone(&self, _seed: u64) -> ArcSampler {
        Arc::new(Self::new(self.data.samples_per_pixel, self.sample_bounds))
    }

    /// This should be called when the rendering algorithm is ready to start
    /// working on a given pixel.
    ///
    /// * `p` - The pixel.
    fn start_pixel(&mut self, p: &Point2i) {
        self.data.start_pixel(p);

        self.gdata.dimension = 0;
        self.gdata.interval_sample_index = self.get_index_for_sample(0);

        // Compute the `array_end_dim` used for aray samples.
        self.gdata.array_end_dim = self.gdata.array_start_dim
            + self.data.sample_array_1d.len() as u16
            + 2 * self.data.sample_array_2d.len() as u16;

        // Compute 1D array samples for `GlobalSampler`.
        let len_1d_sizes = self.data.samples_1d_array_sizes.len();
        for i in 0..len_1d_sizes {
            let n_samples = self.data.samples_1d_array_sizes[i] * self.data.samples_per_pixel;
            for j in 0..n_samples {
                let index = self.get_index_for_sample(j);
                self.data.sample_array_1d[i][j] =
                    self.sample_dimension(index, self.gdata.array_start_dim + i as u16);
            }
        }

        // Compute 2-d array samples for `GlobalSampler`.
        let mut dim = self.gdata.array_start_dim + self.data.samples_1d_array_sizes.len() as u16;
        let len_2d_sizes = self.data.samples_2d_array_sizes.len();
        for i in 0..len_2d_sizes {
            let n_samples = self.data.samples_2d_array_sizes[i] * self.data.samples_per_pixel;
            for j in 0..n_samples {
                let index = self.get_index_for_sample(j);
                self.data.sample_array_2d[i][j] = Point2f::new(
                    self.sample_dimension(index, dim),
                    self.sample_dimension(index, dim + 1),
                );
            }
            dim += 2;
        }

        assert!(self.gdata.array_end_dim == dim);
    }

    /// Returns the sample value for the next dimension of the current sample
    /// vector.
    fn get_1d(&mut self) -> Float {
        if self.gdata.dimension >= self.gdata.array_start_dim
            && self.gdata.dimension < self.gdata.array_end_dim
        {
            self.gdata.dimension = self.gdata.array_end_dim;
        }

        let p = self.sample_dimension(
            self.gdata.interval_sample_index as u64,
            self.gdata.dimension as u16,
        );
        self.gdata.dimension += 1;
        p
    }

    /// Returns the sample value for the next two dimensions of the current
    /// sample vector.
    fn get_2d(&mut self) -> Point2f {
        if self.gdata.dimension + 1 >= self.gdata.array_start_dim
            && self.gdata.dimension < self.gdata.array_end_dim
        {
            self.gdata.dimension = self.gdata.array_end_dim;
        }

        let p = Point2f::new(
            self.sample_dimension(
                self.gdata.interval_sample_index as u64,
                self.gdata.dimension as u16,
            ),
            self.sample_dimension(
                self.gdata.interval_sample_index as u64,
                (self.gdata.dimension + 1) as u16,
            ),
        );
        self.gdata.dimension += 2;
        p
    }

    /// Reset the current sample dimension counter. Returns `true` if
    /// `current_pixel_sample_index` < `samples_per_pixel`; otherwise `false`.
    fn start_next_sample(&mut self) -> bool {
        self.gdata.dimension = 0;
        self.gdata.interval_sample_index =
            self.get_index_for_sample(self.data.current_pixel_sample_index + 1);
        self.data.start_next_sample()
    }

    /// Set the index of the sample in the current pixel to generate next.
    /// Returns `true` if `current_pixel_sample_index` < `samples_per_pixel`;
    /// otherwise `false`.
    ///
    /// * `sample_num` - The sample number.
    fn set_sample_number(&mut self, sample_num: usize) -> bool {
        self.gdata.dimension = 0;
        self.gdata.interval_sample_index = self.get_index_for_sample(sample_num);
        self.data.set_sample_number(sample_num)
    }
}

impl From<(&mut ParamSet, Bounds2i)> for SobolSampler {
    /// Create a `SobolSampler` from given parameter set and sample bounds.
    ///
    /// * `p` - A tuple containing parameter set and sample bounds.
    fn from(p: (&mut ParamSet, Bounds2i)) -> Self {
        let (params, sample_bounds) = p;

        let mut samples_per_pixel = params.find_one_int("pixelsamples", 16) as usize;
        if OPTIONS.quick_render {
            samples_per_pixel = 1;
        }

        Self::new(samples_per_pixel, sample_bounds)
    }
}
