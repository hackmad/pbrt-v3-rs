//! MLT Sampler

use core::geometry::*;
use core::pbrt::*;
use core::rng::*;
use core::sampler::*;

/// MLTSampler maintains three separate sample vectors.
pub const N_SAMPLE_STREAMS: usize = 3;
/// camera subpath
pub const CAMERA_STREAM_INDEX: usize = 0;
/// light subpath
pub const LIGHT_STREAM_INDEX: usize = 1;
/// connection step
pub const CONNECTION_STREAM_INDEX: usize = 2;

/// Implements a sampler responsible for managing primary sample space state vectors, mutations, and acceptance and
/// rejection steps.
pub struct MLTSampler {
    /// The sampler data.
    data: SamplerData,

    /// The random number generator.
    rng: RNG,

    /// Controls the size of “small step” mutations.
    sigma: Float,

    /// Probability of taking a “large step” mutation.
    large_step_probability: Float,

    /// Number of sample streams to request.
    stream_count: usize,

    /// Stores the current sample vector `X`.
    x: Vec<PrimarySample>,

    /// Keeps track of the current Metropolis-Hastings iteration index.
    /// NOTE: Note that iterations with rejected proposals will be excluded from this count.
    current_iteration: usize,

    /// Indicates if a "large step" mutation is taken.
    large_step: bool,

    /// Index of the last iteration where a successful large step took place.
    last_large_step_iteration: usize,

    /// Stream index.
    stream_index: usize,

    /// Index of current sample in the stream.
    sample_index: usize,
}

impl MLTSampler {
    /// Create a new `MLTSampler`.
    ///
    /// * `mutations_per_pixel`    - Number of iterations that MLT (on average!) spends in each pixel.
    /// * `rng_sequence_index`     - Supplies a unique stream index to the internal random number generator.
    /// * `sigma`                  - Controls the size of “small step” mutations.
    /// * `large_step_probability` - Probability of taking a “large step” mutation.
    /// * `stream_count`           - Number of sample streams to request.
    pub fn new(
        mutations_per_pixel: usize,
        rng_sequence_index: u64,
        sigma: Float,
        large_step_probability: Float,
        stream_count: usize,
    ) -> Self {
        Self {
            data: SamplerData::new(mutations_per_pixel),
            rng: RNG::new(rng_sequence_index),
            sigma,
            large_step_probability,
            stream_count,
            x: Vec::new(),
            current_iteration: 0,
            large_step: true,
            last_large_step_iteration: 0,
            stream_index: 0,
            sample_index: 0,
        }
    }

    /// Expands `self.x` as needed and ensures that its contents are in a consistent state.
    ///
    /// * `index` - Index into `self.x`.
    fn ensure_ready(&mut self, index: usize) {
        // Enlarge `MLTSampler::x` if necessary and get current `x_i`.
        if index >= self.x.len() {
            self.x.resize_with(index + 1, Default::default);
        }

        let xi = &mut self.x[index];

        // Reset `x_i` if a large step took place in the meantime.
        if xi.last_modification_iteration < self.last_large_step_iteration {
            xi.value = self.rng.uniform_float();
            xi.last_modification_iteration = self.last_large_step_iteration;
        }

        // Apply remaining sequence of mutations to `sample`.
        xi.backup();
        if self.large_step {
            xi.value = self.rng.uniform_float();
        } else {
            let n_small = self.current_iteration - xi.last_modification_iteration;

            // Apply `n_small` small step mutations.

            // Sample the standard normal distribution `N(0, 1)`.
            let normal_sample = SQRT2 * erf_inv(2.0 * self.rng.uniform_float() - 1.0);

            // Compute the effective standard deviation and apply perturbation to `x_i`.
            let eff_sigma = self.sigma * (n_small as Float).sqrt();
            xi.value += normal_sample * eff_sigma;
            xi.value -= xi.value.floor();
        }
        xi.last_modification_iteration = self.current_iteration;
    }

    /// Call at the beginning of each Metropolis-Hastings iteration. It increases the currentIteration counter and
    /// determines which type of mutation (small or large) should be applied to the sample vector in the current
    /// iteration.
    pub fn start_iteration(&mut self) {
        self.current_iteration += 1;
        self.large_step = self.rng.uniform_float() < self.large_step_probability;
    }

    /// Accept proposed mutation.
    pub fn accept(&mut self) {
        if self.large_step {
            self.last_large_step_iteration = self.current_iteration;
        }
    }

    /// Reject proposed mutation. Restores all `PrimarySamples` modified in the current iteration and reverts the
    /// iteration counter.
    pub fn reject(&mut self) {
        for xi in self.x.iter_mut() {
            if xi.last_modification_iteration == self.current_iteration {
                xi.restore();
            }
        }
        self.current_iteration -= 1;
    }

    /// Indicates that subsequent samples should come from the stream with the given index. It also resets sampleIndex,
    /// the index of the current sample in the stream.
    ///
    /// * `index` - The stream index to start.
    pub fn start_stream(&mut self, index: usize) {
        assert!(index < self.stream_count);
        self.stream_index = index;
        self.sample_index = 0;
    }

    /// Performs corresponding steps through the primary sample vector components. It interleaves the streams into the
    /// global sample vector—in other words, the first `stream_count` dimensions in `x` are respectively used for the
    /// first dimension of each of the streams, and so forth.
    pub fn get_next_index(&mut self) -> usize {
        let next_index = self.stream_index + self.stream_count * self.sample_index;
        self.sample_index += 1;
        next_index
    }
}

impl Sampler for MLTSampler {
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
    /// * `seed` - The seed for the random number generator (ignored).
    fn clone_sampler(&self, seed: u64) -> Box<dyn Sampler> {
        Box::new(Self::new(
            self.data.samples_per_pixel,
            seed,
            self.sigma,
            self.large_step_probability,
            self.stream_count,
        ))
    }

    /// Returns the sample value for the next dimension of the current sample vector.
    fn get_1d(&mut self) -> Float {
        let index = self.get_next_index();
        self.ensure_ready(index);
        self.x[index].value
    }

    /// Returns the sample value for the next two dimensions of the current sample vector.
    fn get_2d(&mut self) -> core::geometry::Point2f {
        Point2f::new(self.get_1d(), self.get_1d())
    }
}

/// Records the current value of a single component of `X` on the interval [0,1).
#[derive(Default)]
pub struct PrimarySample {
    /// Sample value.
    pub value: Float,
    /// Last modification iteration.
    pub last_modification_iteration: usize,
    /// Backed up sample value.
    pub value_backup: Float,
    /// Backed up last modification iteration.
    pub modify_backup: usize,
}

impl PrimarySample {
    /// Backup the sample value.
    fn backup(&mut self) {
        self.value_backup = self.value;
        self.modify_backup = self.last_modification_iteration;
    }

    /// Restore the sample value.
    fn restore(&mut self) {
        self.value = self.value_backup;
        self.last_modification_iteration = self.modify_backup;
    }
}
