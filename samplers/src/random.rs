//! Random Sampler.

use core::geometry::*;
use core::paramset::*;
use core::pbrt::*;
use core::rng::*;
use core::sampler::*;

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
    /// Returns a shared reference underlying `SamplerData`.
    fn get_data(&self) -> &SamplerData {
        &self.data
    }

    /// Returns a mutable reference to underlying `SamplerData`.
    fn get_data_mut(&mut self) -> &mut SamplerData {
        &mut self.data
    }

    /// Generates a new instance of an initial `Sampler` for use by a rendering thread.
    ///
    /// * `seed` - The seed for the random number generator (if any).
    fn clone_sampler(&self, seed: u64) -> Box<dyn Sampler> {
        Box::new(Self::new(self.data.samples_per_pixel, Some(seed)))
    }

    /// This should be called when the rendering algorithm is ready to start working on a given pixel.
    ///
    /// * `p` - The pixel.
    fn start_pixel(&mut self, p: &Point2i) {
        let n = self.data.sample_array_1d.len();
        for i in 0..n {
            for j in 0..self.data.sample_array_1d[i].len() {
                self.data.sample_array_1d[i][j] = self.rng.uniform_float();
            }
        }

        let n = self.data.sample_array_2d.len();
        for i in 0..n {
            for j in 0..self.data.sample_array_2d[i].len() {
                self.data.sample_array_2d[i][j] = Point2f::new(self.rng.uniform_float(), self.rng.uniform_float());
            }
        }

        self.get_data_mut().start_pixel(p);
    }

    /// Returns the sample value for the next dimension of the current sample vector.
    fn get_1d(&mut self) -> Float {
        self.rng.uniform_float()
    }

    /// Returns the sample value for the next two dimensions of the current sample vector.
    fn get_2d(&mut self) -> Point2f {
        Point2f::new(self.rng.uniform_float(), self.rng.uniform_float())
    }
}

impl From<(&ParamSet, Bounds2i)> for RandomSampler {
    /// Create a `RandomSampler` from given parameter set and sample bounds.
    ///
    /// * `p` - A tuple containing parameter set and sample bounds.
    fn from(p: (&ParamSet, Bounds2i)) -> Self {
        let (params, _sample_bounds) = p;
        let samples_per_pixel = params.find_one_int("pixelsamples", 4) as usize;
        Self::new(samples_per_pixel, None)
    }
}
