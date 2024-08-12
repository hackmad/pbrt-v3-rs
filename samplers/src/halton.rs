//! Halton Sampler.

use core::app::OPTIONS;
use core::geometry::*;
use core::low_discrepency::*;
use core::paramset::*;
use core::pbrt::*;
use core::rng::*;
use core::sampler::*;
use std::sync::LazyLock;

/// Maximum resolution for sampling first 2 dimensions.
const K_MAX_RESOLUTION: Int = 128;

/// Precomputed radical inverse permutations.
static RADICAL_INVERSE_PERMUTATIONS: LazyLock<Vec<u16>> = LazyLock::new(|| {
    let mut rng = RNG::default();
    compute_radical_inverse_permutations(&mut rng)
});

/// Implements a low-discrepency sampler using Halton sequences.
pub struct HaltonSampler {
    /// The sampler data.
    data: SamplerData,

    /// The global sampler data.
    gdata: GlobalSamplerData,

    /// Sample bounds.
    sample_bounds: Bounds2i,

    /// The scale factor, either `2^j` or `3^k` for corresponding exponents `j` and `k` stored in `base_exponents`.
    base_scales: Point2<u64>,

    /// The exponents `j` and `k` used to compute the scale factors in `base_scales`.
    base_exponents: Point2<u64>,

    /// Stores the product `2^j * 3^k`. Any particular pixel in the range (0, 0) -> (2^j - 1, 3^k - 1) will be visited
    /// once per `sample_stride`.
    sample_stride: u64,

    /// Multiplicative inverses for `base_scales`.
    mult_inverse: [i64; 2],

    /// Indicates if pixel is sampled at the center.
    sample_at_pixel_center: bool,

    /// Pixel for the current offset.
    pixel_for_offset: Point2i,

    /// Sample index for the first Halton sample for `data.sampler.current_pixel`.
    offset_for_current_pixel: usize,
}

impl HaltonSampler {
    /// Create a new `HaltonSampler`.
    ///
    /// * `samples_per_pixel` - Number of samples per pixel.
    /// * `sample_bounds`     - Sample bounds.
    /// * `sample_at_center`  - Indicates whether or not to jitter each sample's center point (default to false).
    pub fn new(samples_per_pixel: usize, sample_bounds: Bounds2i, sample_at_center: bool) -> Self {
        // Find radical inverse, base scales and exponents that cover sampling area.
        let res = sample_bounds.p_max - sample_bounds.p_min;
        let mut base_scales = Point2::<u64>::default();
        let mut base_exponents = Point2::<u64>::default();

        for i in 0..2 {
            let base = if i == 0 { 2_u64 } else { 3_u64 };
            let mut scale = 1_u64;
            let mut exp = 0_u64;
            while (scale as Int) < min(res[i], K_MAX_RESOLUTION) {
                scale *= base;
                exp += 1;
            }
            base_scales[i] = scale;
            base_exponents[i] = exp;
        }

        // Compute stride in samples for visiting each pixel area.
        let sample_stride = base_scales[0] * base_scales[1];

        // Compute multiplicative inverses for `base_scales`.
        let mult_inverse = [
            multiplicative_inverse(base_scales[1] as i64, base_scales[0] as i64) as i64,
            multiplicative_inverse(base_scales[0] as i64, base_scales[1] as i64) as i64,
        ];

        Self {
            data: SamplerData::new(samples_per_pixel),
            gdata: GlobalSamplerData::new(),
            sample_at_pixel_center: sample_at_center,
            pixel_for_offset: Point2i::new(Int::MAX, Int::MAX),
            offset_for_current_pixel: 0,
            sample_bounds,
            base_scales,
            base_exponents,
            sample_stride,
            mult_inverse,
        }
    }

    /// Returns the radical inverse permutation for a given dimension.
    ///
    /// * `dim` - Dimension.
    fn permutation_for_dimension(&self, dim: u16) -> &[u16] {
        assert!(
            (dim as usize) <= PRIME_TABLE_SIZE,
            "HaltonSampler can only sample {} dimensions",
            PRIME_TABLE_SIZE
        );
        &RADICAL_INVERSE_PERMUTATIONS[PRIME_SUMS[dim as usize]..]
    }

    /// Performs the inverse mapping from the current pixel and given sample index to a global index into the overall
    /// set of sample vectors.
    ///
    /// * `sample_num` - The sample number.
    fn get_index_for_sample(&mut self, sample_num: usize) -> u64 {
        if self.data.current_pixel != self.pixel_for_offset {
            // Compute Halton sample offset for _currentPixel_
            self.offset_for_current_pixel = 0;
            if self.sample_stride > 1 {
                let pm = Point2i::new(
                    rem(self.data.current_pixel[0], K_MAX_RESOLUTION),
                    rem(self.data.current_pixel[1], K_MAX_RESOLUTION),
                );
                for i in 0..2 {
                    let dim_offset = if i == 0 {
                        inverse_radical_inverse(2, pm[i] as u64, self.base_exponents[i])
                    } else {
                        inverse_radical_inverse(3, pm[i] as u64, self.base_exponents[i])
                    };
                    let offset = dim_offset * (self.sample_stride / self.base_scales[i]) * self.mult_inverse[i] as u64;
                    self.offset_for_current_pixel += offset as usize;
                }

                self.offset_for_current_pixel %= self.sample_stride as usize;
            }

            self.pixel_for_offset = self.data.current_pixel;
        }

        (self.offset_for_current_pixel + sample_num * self.sample_stride as usize) as u64
    }

    /// Returns the sample value for the given dimension of the index^th sample vector in the sequence.
    ///
    /// * `index` - Index of the sample.
    /// * `dim`   - Dimension.
    fn sample_dimension(&mut self, index: u64, dim: u16) -> Float {
        if self.sample_at_pixel_center && (dim == 0 || dim == 1) {
            0.5
        } else if dim == 0 {
            radical_inverse(dim, index >> self.base_exponents[0])
        } else if dim == 1 {
            radical_inverse(dim, index / self.base_scales[1])
        } else {
            scrambled_radical_inverse(dim, index, self.permutation_for_dimension(dim))
        }
    }
}

impl Sampler for HaltonSampler {
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
    fn clone_sampler(&self, _seed: u64) -> Box<dyn Sampler> {
        Box::new(Self::new(
            self.data.samples_per_pixel,
            self.sample_bounds,
            self.sample_at_pixel_center,
        ))
    }

    /// This should be called when the rendering algorithm is ready to start working on a given pixel.
    ///
    /// * `p` - The pixel.
    fn start_pixel(&mut self, p: &Point2i) {
        self.data.start_pixel(p);

        self.gdata.dimension = 0;
        self.gdata.interval_sample_index = self.get_index_for_sample(0);

        // Compute the `array_end_dim` used for array samples.
        self.gdata.array_end_dim = self.gdata.array_end_dim
            + self.data.sample_array_1d.len() as u16
            + 2 * self.data.sample_array_2d.len() as u16;

        // Compute 1D array samples for `GlobalSampler`.
        let len_1d_sizes = self.data.samples_1d_array_sizes.len();
        for i in 0..len_1d_sizes {
            let n_samples = self.data.samples_1d_array_sizes[i] * self.data.samples_per_pixel;
            for j in 0..n_samples {
                let index = self.get_index_for_sample(j);
                self.data.sample_array_1d[i][j] = self.sample_dimension(index, self.gdata.array_start_dim + i as u16);
            }
        }

        // Compute 2-d array samples for `GlobalSampler`.
        let mut dim = self.gdata.array_start_dim + self.data.samples_1d_array_sizes.len() as u16;
        let len_2d_sizes = self.data.samples_2d_array_sizes.len();
        for i in 0..len_2d_sizes {
            let n_samples = self.data.samples_2d_array_sizes[i] * self.data.samples_per_pixel;
            for j in 0..n_samples {
                let index = self.get_index_for_sample(j);
                self.data.sample_array_2d[i][j] =
                    Point2f::new(self.sample_dimension(index, dim), self.sample_dimension(index, dim + 1));
            }
            dim += 2;
        }

        assert!(self.gdata.array_end_dim == dim);
    }

    /// Returns the sample value for the next dimension of the current sample vector.
    fn get_1d(&mut self) -> Float {
        if self.gdata.dimension >= self.gdata.array_start_dim && self.gdata.dimension < self.gdata.array_end_dim {
            self.gdata.dimension = self.gdata.array_end_dim;
        }

        let p = self.sample_dimension(self.gdata.interval_sample_index as u64, self.gdata.dimension as u16);
        self.gdata.dimension += 1;
        p
    }

    /// Returns the sample value for the next two dimensions of the current sample vector.
    fn get_2d(&mut self) -> Point2f {
        if self.gdata.dimension + 1 >= self.gdata.array_start_dim && self.gdata.dimension < self.gdata.array_end_dim {
            self.gdata.dimension = self.gdata.array_end_dim;
        }

        let p = Point2f::new(
            self.sample_dimension(self.gdata.interval_sample_index as u64, self.gdata.dimension as u16),
            self.sample_dimension(
                self.gdata.interval_sample_index as u64,
                (self.gdata.dimension + 1) as u16,
            ),
        );
        self.gdata.dimension += 2;
        p
    }

    /// Reset the current sample dimension counter. Returns `true` if `current_pixel_sample_index` <
    /// `samples_per_pixel`; otherwise `false`.
    fn start_next_sample(&mut self) -> bool {
        self.gdata.dimension = 0;
        self.gdata.interval_sample_index = self.get_index_for_sample(self.data.current_pixel_sample_index + 1);
        self.data.start_next_sample()
    }

    /// Set the index of the sample in the current pixel to generate next. Returns `true` if
    /// `current_pixel_sample_index` < `samples_per_pixel`; otherwise `false`.
    ///
    /// * `sample_num` - The sample number.
    fn set_sample_number(&mut self, sample_num: usize) -> bool {
        self.gdata.dimension = 0;
        self.gdata.interval_sample_index = self.get_index_for_sample(sample_num);
        self.data.set_sample_number(sample_num)
    }
}

impl From<(&ParamSet, Bounds2i)> for HaltonSampler {
    /// Create a `HaltonSampler` from given parameter set and sample bounds.
    ///
    /// * `p` - A tuple containing parameter set and sample bounds.
    fn from(p: (&ParamSet, Bounds2i)) -> Self {
        let (params, sample_bounds) = p;

        let mut samples_per_pixel = params.find_one_int("pixelsamples", 16) as usize;
        if OPTIONS.quick_render {
            samples_per_pixel = 1;
        }

        let sample_at_center = params.find_one_bool("samplepixelcenter", false);

        Self::new(samples_per_pixel, sample_bounds, sample_at_center)
    }
}

/// Returns the GCD of two numbers.
///
/// * `a` - First number.
/// * `b` - Second number.
fn extended_gcd(a: u64, b: u64) -> (i64, i64) {
    if b == 0 {
        (1, 0)
    } else {
        let d = (a / b) as i64;
        let (xp, yp) = extended_gcd(b, a % b);
        (yp, xp - (d * yp))
    }
}

/// Calculate the multiplicative inverse `b` of `a` with respect to modulus `n` such that `(a * b) mod n = 1`.
///
/// * `a` - Number.
/// * `n` - Modulus.
fn multiplicative_inverse(a: i64, n: i64) -> u64 {
    let (x, _y) = extended_gcd(a as u64, n as u64);
    rem(x, n) as u64
}
