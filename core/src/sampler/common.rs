//! Common

use super::*;

/// Stores the sampler data and implements common functionality for all samplers.
#[derive(Clone, Default)]
pub struct SamplerData {
    /// Number of samples generated for each pixel.
    pub samples_per_pixel: usize,

    /// Coordinates of current pixel being generated.
    pub current_pixel: Point2i,

    /// Sample number of the pixel currently being generated.
    pub current_pixel_sample_index: usize,

    /// Stores sizes of requested 1D sample arrays.
    pub samples_1d_array_sizes: Vec<usize>,

    /// Stores sizes of requested 2D sample arrays.
    pub samples_2d_array_sizes: Vec<usize>,

    /// Stores the `n` requested 1D samples for each pixel.
    pub sample_array_1d: Vec<Vec<Float>>,

    /// Stores the `n` requested 2D samples for each pixel.
    pub sample_array_2d: Vec<Vec<Point2f>>,

    /// Tracks index of the next element in 1D array. This is reset to 0 when
    /// a new samplpixel starts or the sample number in current pixel changes.
    pub array_1d_offset: usize,

    /// Tracks index of the next element in 2D array. This is reset to 0 when
    /// a new samplpixel starts or the sample number in current pixel changes.
    pub array_2d_offset: usize,
}

impl SamplerData {
    /// Create a new `SamplerData` instance.
    ///
    /// * `samples_per_pixel` - Number of samples to generate for each pixel.
    pub fn new(samples_per_pixel: usize) -> Self {
        Self {
            samples_per_pixel,
            current_pixel: Point2i::default(),
            current_pixel_sample_index: 0,
            samples_1d_array_sizes: vec![],
            samples_2d_array_sizes: vec![],
            sample_array_1d: vec![],
            sample_array_2d: vec![],
            array_1d_offset: 0,
            array_2d_offset: 0,
        }
    }

    /// This should be called when the rendering algorithm is ready to start
    /// working on a given pixel.
    ///
    /// * `p` - The pixel.
    pub fn start_pixel(&mut self, p: &Point2i) {
        self.current_pixel = *p;
        self.current_pixel_sample_index = 0;

        // Reset array offsets for next pixel sample.
        self.array_1d_offset = 0;
        self.array_2d_offset = 0;
    }

    /// This should be called before rendering begins when an array of 1D
    /// samples is required.
    ///
    /// * `n` - The number of samples.
    pub fn request_1d_array(&mut self, n: usize) {
        self.samples_1d_array_sizes.push(n);
        self.sample_array_1d
            .push(Vec::<Float>::with_capacity(n * self.samples_per_pixel));
    }

    /// This should be called before rendering begins when an array of 1D
    /// samples is required.
    ///
    /// * `n` - The number of samples.
    pub fn request_2d_array(&mut self, n: usize) {
        self.samples_2d_array_sizes.push(n);
        self.sample_array_2d
            .push(Vec::<Point2f>::with_capacity(n * self.samples_per_pixel));
    }

    /// Get an array of 1D samples.
    ///
    /// * `n` - The number of samples.
    pub fn get_1d_array(&mut self, n: usize) -> Vec<Float> {
        if self.array_1d_offset == self.sample_array_1d.len() {
            vec![]
        } else {
            assert!(self.samples_1d_array_sizes[self.array_1d_offset] == n);
            assert!(self.current_pixel_sample_index < self.samples_per_pixel);

            let array = &self.sample_array_1d[self.array_1d_offset];
            self.array_1d_offset += 1;

            let i = self.current_pixel_sample_index * n;
            let m = i + n;
            array[i..m].to_vec()
        }
    }

    /// Get an array of 2D samples.
    ///
    /// * `n` - The number of samples.
    pub fn get_2d_array(&mut self, n: usize) -> Vec<Point2f> {
        if self.array_2d_offset == self.sample_array_2d.len() {
            vec![]
        } else {
            assert!(self.samples_2d_array_sizes[self.array_2d_offset] == n);
            assert!(self.current_pixel_sample_index < self.samples_per_pixel);

            let array = &self.sample_array_2d[self.array_2d_offset];
            self.array_2d_offset += 1;

            let i = self.current_pixel_sample_index * n;
            let m = i + n;
            array[i..m].to_vec()
        }
    }

    /// Reset the current sample dimension counter. Returns `true` if
    /// `current_pixel_sample_index` < `samples_per_pixel`; otherwise `false`.
    pub fn start_next_sample(&mut self) -> bool {
        // Reset array offsets for next pixel sample.
        self.array_1d_offset = 0;
        self.array_2d_offset = 0;
        self.current_pixel_sample_index += 1;
        self.current_pixel_sample_index < self.samples_per_pixel
    }

    /// Set the index of the sample in the current pixel to generate next.
    /// Returns `true` if `current_pixel_sample_index` < `samples_per_pixel`;
    /// otherwise `false`.
    ///
    /// * `sample_num` - The sample number.
    pub fn set_sample_number(&mut self, sample_num: usize) -> bool {
        // Reset array offsets for next pixel sample.
        self.array_1d_offset = 0;
        self.array_2d_offset = 0;
        self.current_pixel_sample_index = sample_num;
        self.current_pixel_sample_index < self.samples_per_pixel
    }

    /// Returns the index of the samle in the current pixel.
    pub fn current_sample_number(&self) -> usize {
        self.current_pixel_sample_index
    }
}

/// Stores the data for samplers that are not pixel based and can generate
/// samples that are spread across the entire image. This can support samplers
/// that can generate samples for image tiles.
#[derive(Clone)]
pub struct GlobalSamplerData {
    /// Sample dimension.
    pub dimension: u16,

    /// The index of the sample in the current pixel.
    pub interval_sample_index: u64,

    /// The first dimensions up to this are devoted to regular 1D and 2D
    /// samples. Subsequent dimensions are devoted to first 1D, then 2D
    /// array samples up to and not included `array_end_dim`.
    pub array_start_dim: u16,

    /// Higher dimension samples start here and are used for non-array 1D and
    /// 2D samples.
    pub array_end_dim: u16,
}

impl GlobalSamplerData {
    /// Create a new `GlobalSamplerData` instance.
    ///
    /// * `samples_per_pixel` - Number of samples to generate for each pixel.
    pub fn new() -> Self {
        Self {
            dimension: 0,
            interval_sample_index: 0,
            array_start_dim: 5,
            array_end_dim: 5,
        }
    }
}

impl Default for GlobalSamplerData {
    /// Returns the "default value" for `GlobalSamplerData`.
    fn default() -> Self {
        Self::new()
    }
}
