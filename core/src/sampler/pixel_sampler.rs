//! Pixel Sampler.

use super::*;
use std::sync::Arc;

/// Implementation for generating all sample values for all sample vectors of
/// a pixel at a time.
///
/// If you use this in your own sampler which implements `Sampler` interface,
/// make sure to call functions in this implementation if you do not wish to
/// implement them youself. Since this implements `Sampler` interface, it will
/// forward the call up to default implementation if its not implemented here
/// preserving the inheritence structure of the original PBRT implementation.
#[derive(Clone)]
pub struct PixelSampler {
    /// The common sampler data.
    pub data: SamplerData,

    /// Vector of 1D sample values.
    pub samples_1d: Vec<Vec<Float>>,

    /// Vector of 2D sample values.
    pub samples_2d: Vec<Vec<Point2f>>,

    /// Offset into the `samples_1d` array for the current pixel sample. It must
    /// be reset to 0 at the start of each new sample.
    pub current_1d_dimension: usize,

    /// Offset into the `samples_2d` array for the current pixel sample. It must
    /// be reset to 0 at the start of each new sample.
    pub current_2d_dimension: usize,

    /// The random number generator.
    pub rng: RNG,
}

impl PixelSampler {
    /// Create a new `PixelSampler`.
    ///
    /// * `samples_per_pixel`    - Number of samples per pixel.
    /// * `n_sampled_dimensions` - Number of dimensions for sampling.
    /// * `seed`                 - Optional seed for the random number generator.
    pub fn new(samples_per_pixel: usize, n_sampled_dimensions: usize, seed: Option<u64>) -> Self {
        let mut samples_1d = Vec::<Vec<Float>>::with_capacity(n_sampled_dimensions);
        let mut samples_2d = Vec::<Vec<Point2f>>::with_capacity(n_sampled_dimensions);

        for _i in 0..n_sampled_dimensions {
            samples_1d.push(vec![0.0; samples_per_pixel]);
            samples_2d.push(vec![Point2f::ZERO; samples_per_pixel]);
        }

        let rng = match seed {
            Some(s) => RNG::new(s),
            None => RNG::default(),
        };

        Self {
            data: SamplerData::new(samples_per_pixel),
            samples_1d,
            samples_2d,
            current_1d_dimension: 0,
            current_2d_dimension: 0,
            rng,
        }
    }
}

impl Sampler for PixelSampler {
    /// Returns the underlying `SamplerData`.
    fn get_data(&mut self) -> &mut SamplerData {
        &mut self.data
    }

    /// Generates a new instance of an initial `Sampler` for use by a rendering
    /// thread.
    ///
    /// * `seed` - The seed for the random number generator (if any).
    fn clone(&self, seed: u64) -> ArcSampler {
        Arc::new(Self::new(
            self.data.samples_per_pixel,
            self.samples_1d.len(),
            Some(seed),
        ))
    }

    /// Returns the sample value for the next dimension of the current sample
    /// vector.
    fn get_1d(&mut self) -> Float {
        assert!(self.data.current_pixel_sample_index < self.data.samples_per_pixel);
        if self.current_1d_dimension < self.samples_1d.len() {
            let r =
                self.samples_1d[self.current_1d_dimension][self.data.current_pixel_sample_index];
            self.current_1d_dimension += 1;
            r
        } else {
            self.rng.uniform_float()
        }
    }

    /// Returns the sample value for the next two dimensions of the current
    /// sample vector.
    fn get_2d(&mut self) -> Point2f {
        assert!(self.data.current_pixel_sample_index < self.data.samples_per_pixel);
        if self.current_2d_dimension < self.samples_2d.len() {
            let r =
                self.samples_2d[self.current_2d_dimension][self.data.current_pixel_sample_index];
            self.current_2d_dimension += 1;
            r
        } else {
            Point2f::new(self.rng.uniform_float(), self.rng.uniform_float())
        }
    }

    /// Reset the current sample dimension counter. Returns `true` if
    /// `current_pixel_sample_index` < `samples_per_pixel`; otherwise `false`.
    fn start_next_sample(&mut self) -> bool {
        self.current_1d_dimension = 0;
        self.current_2d_dimension = 0;
        self.data.start_next_sample()
    }

    /// Set the index of the sample in the current pixel to generate next.
    /// Returns `true` if `current_pixel_sample_index` < `samples_per_pixel`;
    /// otherwise `false`.
    ///
    /// * `sample_num` - The sample number.
    fn set_sample_number(&mut self, sample_num: usize) -> bool {
        self.current_1d_dimension = 0;
        self.current_2d_dimension = 0;
        self.data.set_sample_number(sample_num)
    }
}
