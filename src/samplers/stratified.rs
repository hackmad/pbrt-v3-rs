//! Stratified Sampler.

use super::SamplerProps;
use crate::core::app::OPTIONS;
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::sampler::*;
use crate::core::sampling::*;
use std::sync::Arc;

// Implements a stratified sampler that subdivides pixel areas into non-overlapping
/// rectangular regions, called strata, and geenrates a single sample inside each
/// region.
pub struct StratifiedSampler {
    /// Use a pixel sampler.
    sampler: PixelSampler,

    /// Number of samples in horizontal direction.
    x_pixel_samples: usize,

    /// Number of samples in vertical direction.
    y_pixel_samples: usize,

    /// Indicates whether or not to jitter each sample's center point.
    /// In unjittered mode is is uniform sampling and is not useful for
    /// high quality rendering but mostly for comparing sampling techniques.
    jitter_samples: bool,
}

impl StratifiedSampler {
    /// Create a new `StratifiedSampler`.
    ///
    /// * `x_pixel_samples`      - Number of samples in horizontal direction.
    /// * `y_pixel_samples`      - Number of samples in vertical direction.
    /// * `jitter_samples`       - Indicates whether or not to jitter each
    ///                            sample's center point.
    /// * `n_sampled_dimensions` - Number of dimensions for sampling.
    /// * `seed`                 - Optional seed for the random number generator.
    pub fn new(
        x_pixel_samples: usize,
        y_pixel_samples: usize,
        jitter_samples: bool,
        n_sampled_dimensions: usize,
        seed: Option<u64>,
    ) -> Self {
        let samples_per_pixel = x_pixel_samples * y_pixel_samples;
        Self {
            sampler: PixelSampler::new(samples_per_pixel, n_sampled_dimensions, seed),
            x_pixel_samples,
            y_pixel_samples,
            jitter_samples,
        }
    }
}

impl Sampler for StratifiedSampler {
    /// Returns the underlying `SamplerData`.
    fn get_data(&mut self) -> &mut SamplerData {
        &mut self.sampler.data
    }

    /// Generates a new instance of an initial `Sampler` for use by a rendering
    /// thread.
    ///
    /// * `seed` - The seed for the random number generator (if any).
    fn clone(&self, seed: u64) -> ArcSampler {
        Arc::new(Self::new(
            self.x_pixel_samples,
            self.y_pixel_samples,
            self.jitter_samples,
            self.sampler.samples_1d.len(),
            Some(seed),
        ))
    }

    /// This should be called when the rendering algorithm is ready to start
    /// working on a given pixel.
    ///
    /// * `p` - The pixel.
    fn start_pixel(&mut self, p: &Point2i) {
        let samples_per_pixel = self.sampler.data.samples_per_pixel;

        // Generate single stratified samples for the pixel.
        for i in 0..self.sampler.samples_1d.len() {
            let mut samples = stratified_sample_1d(
                &mut self.sampler.rng,
                samples_per_pixel,
                self.jitter_samples,
            );
            self.sampler.rng.shuffle(&mut samples, samples_per_pixel, 1);
            for j in 0..samples_per_pixel {
                self.sampler.samples_1d[i][j] = samples[j];
            }
        }

        for i in 0..self.sampler.samples_2d.len() {
            let mut samples = stratified_sample_2d(
                &mut self.sampler.rng,
                self.x_pixel_samples,
                self.y_pixel_samples,
                self.jitter_samples,
            );
            self.sampler.rng.shuffle(&mut samples, samples_per_pixel, 1);
            for j in 0..samples_per_pixel {
                self.sampler.samples_2d[i][j] = samples[j];
            }
        }

        // Generate arrays of stratified samples for the pixel.
        for i in 0..self.sampler.data.samples_1d_array_sizes.len() {
            for j in 0..samples_per_pixel {
                let count = self.sampler.data.samples_1d_array_sizes[i];
                let mut samples =
                    stratified_sample_1d(&mut self.sampler.rng, count, self.jitter_samples);
                self.sampler.rng.shuffle(&mut samples, count, 1);
                for k in 0..count {
                    self.sampler.data.sample_array_1d[i][j * count + k];
                }
            }
        }

        for i in 0..self.sampler.data.samples_2d_array_sizes.len() {
            for j in 0..samples_per_pixel {
                let count = self.sampler.data.samples_2d_array_sizes[i];
                let samples = latin_hypercube(&mut self.sampler.rng, count, 2);
                for k in 0..count * 2 {
                    self.sampler.data.sample_array_2d[i][j * count + k] =
                        Point2f::new(samples[k], samples[k + 1]);
                }
            }
        }

        self.get_data().start_pixel(p);
    }

    /// Returns the sample value for the next dimension of the current sample
    /// vector.
    fn get_1d(&mut self) -> Float {
        self.sampler.get_1d()
    }

    /// Returns the sample value for the next two dimensions of the current
    /// sample vector.
    fn get_2d(&mut self) -> Point2f {
        self.sampler.get_2d()
    }

    /// Reset the current sample dimension counter. Returns `true` if
    /// `current_pixel_sample_index` < `samples_per_pixel`; otherwise `false`.
    fn start_next_sample(&mut self) -> bool {
        self.sampler.start_next_sample()
    }

    /// Set the index of the sample in the current pixel to generate next.
    /// Returns `true` if `current_pixel_sample_index` < `samples_per_pixel`;
    /// otherwise `false`.
    ///
    /// * `sample_num` - The sample number.
    fn set_sample_number(&mut self, sample_num: usize) -> bool {
        self.sampler.set_sample_number(sample_num)
    }
}

impl From<&mut SamplerProps> for StratifiedSampler {
    /// Create a `StratifiedSampler` from `SamplerProps`.
    ///
    /// * `props` - Sampler creation properties.
    fn from(props: &mut SamplerProps) -> Self {
        let mut x_samples = props.params.find_one_int("xsamples", 4) as usize;
        let mut y_samples = props.params.find_one_int("ysamples", 4) as usize;
        if OPTIONS.quick_render {
            x_samples = 1;
            y_samples = 1;
        }

        let jitter = props.params.find_one_bool("jitter", true);
        let sd = props.params.find_one_int("dimensions", 4) as usize;

        Self::new(x_samples, y_samples, jitter, sd, None)
    }
}
