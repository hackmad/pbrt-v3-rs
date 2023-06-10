//! Maximized Minimal Distance Sampler.

use core::app::options;
use core::geometry::*;
use core::low_discrepency::*;
use core::paramset::*;
use core::pbrt::*;
use core::sampler::*;

/// Implements a sampler that can generate (0, 2)-sequences but also maximizes the distance between samples.
pub struct MaxMinDistSampler {
    /// Use a pixel sampler.
    sampler: PixelSampler,

    /// Stores the specialized generator matrix.
    c_pixel: [u32; 32],
}

impl MaxMinDistSampler {
    /// Create a new `MaxMinDistSampler`.
    ///
    /// * `samples_per_pixel`    - Number of samples per pixel.
    /// * `n_sampled_dimensions` - Number of dimensions for sampling.
    /// * `seed`                 - Optional seed for the random number generator.
    pub fn new(samples_per_pixel: usize, n_sampled_dimensions: usize, seed: Option<u64>) -> Self {
        let mut spp = samples_per_pixel;

        let mut c_index = spp.log2int();
        let max_c_index = C_MAX_MIN_DIST.len();

        if c_index as usize >= max_c_index {
            warn!(
                "No more than {} samples per pixel are supported with \
                MaxMinDistSampler. Rounding down.",
                spp
            );
            spp = (1 << max_c_index) - 1;
        }
        if !spp.is_power_of_two() {
            let spp2 = samples_per_pixel.next_power_of_two();
            warn!(
                "Non power-of-two {} sample count for MaxMinDistSampler. \
                Rounded up to {}.",
                spp, spp2
            );
            spp = spp2;
        }

        c_index = spp.log2int();
        assert!(c_index >= 0_i64 && c_index < max_c_index as i64);

        Self {
            sampler: PixelSampler::new(spp, n_sampled_dimensions, seed),
            c_pixel: C_MAX_MIN_DIST[c_index as usize],
        }
    }
}

impl Sampler for MaxMinDistSampler {
    /// Returns a shared reference underlying `SamplerData`.
    fn get_data(&self) -> &SamplerData {
        &self.sampler.data
    }

    /// Returns a mutable reference to underlying `SamplerData`.
    fn get_data_mut(&mut self) -> &mut SamplerData {
        &mut self.sampler.data
    }

    /// Generates a new instance of an initial `Sampler` for use by a rendering thread.
    ///
    /// * `seed` - The seed for the random number generator (if any).
    fn clone_sampler(&self, seed: u64) -> Box<dyn Sampler> {
        Box::new(Self::new(
            self.sampler.data.samples_per_pixel,
            self.sampler.samples_1d.len(),
            Some(seed),
        ))
    }

    /// This should be called when the rendering algorithm is ready to start working on a given pixel.
    ///
    /// * `p` - The pixel.
    fn start_pixel(&mut self, p: &Point2i) {
        let samples_per_pixel = self.sampler.data.samples_per_pixel;

        let inv_spp = 1.0 / samples_per_pixel as Float;

        for i in 0..samples_per_pixel {
            self.sampler.samples_2d[0][i] = Point2f::new(
                i as Float * inv_spp,
                sample_generator_matrix(&self.c_pixel, i as u32, 0),
            );
        }
        self.sampler
            .rng
            .shuffle(&mut self.sampler.samples_2d[0], samples_per_pixel, 1);

        // Generate remaining samples for `MaxMinDistSampler`.
        let n = self.sampler.samples_1d.len();
        for i in 0..n {
            van_der_corput(
                1,
                samples_per_pixel,
                &mut self.sampler.samples_1d[i],
                &mut self.sampler.rng,
            );
        }

        let n = self.sampler.samples_2d.len();
        for i in 1..n {
            sobol_2d(
                1,
                samples_per_pixel,
                &mut self.sampler.samples_2d[i],
                &mut self.sampler.rng,
            );
        }

        let n = self.sampler.data.samples_1d_array_sizes.len();
        for i in 0..n {
            van_der_corput(
                self.sampler.data.samples_1d_array_sizes[i],
                samples_per_pixel,
                &mut self.sampler.data.sample_array_1d[i],
                &mut self.sampler.rng,
            );
        }

        let n = self.sampler.data.samples_2d_array_sizes.len();
        for i in 0..n {
            sobol_2d(
                self.sampler.data.samples_2d_array_sizes[i],
                samples_per_pixel,
                &mut self.sampler.data.sample_array_2d[i],
                &mut self.sampler.rng,
            );
        }

        self.get_data_mut().start_pixel(p);
    }

    /// Returns the sample value for the next dimension of the current sample vector.
    fn get_1d(&mut self) -> Float {
        self.sampler.get_1d()
    }

    /// Returns the sample value for the next two dimensions of the current sample vector.
    fn get_2d(&mut self) -> Point2f {
        self.sampler.get_2d()
    }

    /// Rounds up to power of 2.
    ///
    /// * `n` - The integer value to round.
    fn round_count(&self, n: usize) -> usize {
        n.next_power_of_two()
    }

    /// Reset the current sample dimension counter. Returns `true` if `current_pixel_sample_index` <
    /// `samples_per_pixel`; otherwise `false`.
    fn start_next_sample(&mut self) -> bool {
        self.sampler.start_next_sample()
    }

    /// Set the index of the sample in the current pixel to generate next. Returns `true` if
    /// `current_pixel_sample_index` < `samples_per_pixel`; otherwise `false`.
    ///
    /// * `sample_num` - The sample number.
    fn set_sample_number(&mut self, sample_num: usize) -> bool {
        self.sampler.set_sample_number(sample_num)
    }
}

impl From<(&ParamSet, Bounds2i)> for MaxMinDistSampler {
    /// Create a `MaxMinDistSampler` from given parameter set and sample bounds.
    ///
    /// * `p` - A tuple containing parameter set and sample bounds.
    fn from(p: (&ParamSet, Bounds2i)) -> Self {
        let (params, _sample_bounds) = p;

        let mut samples_per_pixel = params.find_one_int("pixelsamples", 16) as usize;
        if options().quick_render {
            samples_per_pixel = 1;
        }

        let sd = params.find_one_int("dimensions", 4) as usize;

        Self::new(samples_per_pixel, sd, None)
    }
}
