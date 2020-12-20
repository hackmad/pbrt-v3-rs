//! Random Sampler.

use super::{ArcSampler, Float, Point2f, Point2i, Sampler, SamplerData, UniformRandom, RNG};
use std::sync::Arc;

/// Implements a sampler that uses a PRNG to generate uniformly random samples.
pub struct RandomSampler {
    /// The common sampler data.
    pub data: SamplerData,

    /// The random number generator.
    pub rng: RNG,
}

impl RandomSampler {
    /// Create a new `RandomSampler`.
    ///
    /// * `samples_per_pixel` - Number of samples to generate for each pixel.
    /// * `seed`              - Optional seed for the random number generator.
    pub fn new(samples_per_pixel: usize, seed: Option<u64>) -> Self {
        Self {
            data: SamplerData::new(samples_per_pixel),
            rng: match seed {
                Some(s) => RNG::new(s),
                None => RNG::default(),
            },
        }
    }
}

impl Sampler for RandomSampler {
    /// Returns the underlying `SamplerData`.
    fn get_data(&mut self) -> &mut SamplerData {
        &mut self.data
    }

    /// Generates a new instance of an initial `Sampler` for use by a rendering
    /// thread.
    ///
    /// * `seed` - The seed for the random number generator (if any).
    fn clone(&self, seed: u64) -> ArcSampler {
        Arc::new(Self::new(self.data.samples_per_pixel, Some(seed)))
    }

    /// This should be called when the rendering algorithm is ready to start
    /// working on a given pixel.
    ///
    /// * `p` - The pixel.
    fn start_pixel(&mut self, p: &Point2i) {
        let n = self.data.sample_array_1d.len();
        for i in 0..n {
            for j in 0..self.data.sample_array_1d[i].len() {
                self.data.sample_array_1d[i][j] = self.rng.uniform();
            }
        }

        let n = self.data.sample_array_2d.len();
        for i in 0..n {
            for j in 0..self.data.sample_array_2d[i].len() {
                self.data.sample_array_2d[i][j] =
                    Point2f::new(self.rng.uniform(), self.rng.uniform());
            }
        }

        self.get_data().start_pixel(p);
    }

    /// Returns the sample value for the next dimension of the current sample
    /// vector.
    fn get_1d(&mut self) -> Float {
        self.rng.uniform()
    }

    /// Returns the sample value for the next two dimensions of the current
    /// sample vector.
    fn get_2d(&mut self) -> Point2f {
        Point2f::new(self.rng.uniform(), self.rng.uniform())
    }
}
