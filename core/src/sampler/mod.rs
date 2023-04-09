//! Sampler

mod common;
mod pixel_sampler;

use crate::camera::*;
use crate::geometry::*;
use crate::pbrt::*;
use crate::rng::*;
use std::sync::Arc;

// Re-export
pub use common::*;
pub use pixel_sampler::*;

/// Sampler interface.
pub trait Sampler {
    /// Returns a shared reference underlying `SamplerData`.
    fn get_data(&self) -> &SamplerData;

    /// Returns a mutable reference to underlying `SamplerData`.
    fn get_data_mut(&mut self) -> &mut SamplerData;

    /// Generates a new instance of an initial `Sampler` for use by a rendering thread.
    ///
    /// * `seed` - The seed for the random number generator (if any).
    fn clone_sampler(&self, seed: u64) -> Box<dyn Sampler>;

    /// This should be called when the rendering algorithm is ready to start working on a given pixel.
    ///
    /// * `p` - The pixel.
    fn start_pixel(&mut self, p: &Point2i) {
        self.get_data_mut().start_pixel(p);
    }

    /// Returns the sample value for the next dimension of the current sample vector.
    fn get_1d(&mut self) -> Float;

    /// Returns the sample value for the next two dimensions of the current sample vector.
    fn get_2d(&mut self) -> Point2f;

    /// Returns an initialized `CameraSample` for a given pixel.
    ///
    /// * `p_raster` - The pixel.
    fn get_camera_sample(&mut self, p_raster: &Point2i) -> CameraSample {
        let film_sample = self.get_2d();
        let p_film = Point2f::new(p_raster.x as Float + film_sample.x, p_raster.y as Float + film_sample.y);

        let time = self.get_1d();
        let p_lens = self.get_2d();

        CameraSample::new(p_film, p_lens, time)
    }

    /// This should be called before rendering begins when an array of 1D samples is required.
    ///
    /// * `n` - The number of samples.
    fn request_1d_array(&mut self, n: usize) {
        assert!(self.round_count(n) == n);
        self.get_data_mut().request_1d_array(n);
    }

    /// This should be called before rendering begins when an array of 1D samples is required.
    ///
    /// * `n` - The number of samples.
    fn request_2d_array(&mut self, n: usize) {
        assert!(self.round_count(n) == n);
        self.get_data_mut().request_2d_array(n);
    }

    /// Returns nearest interger based on some criteria (e.g. nearest power of two). The default implementation simply
    /// returns the given value.
    ///
    /// * `n` - The integer value to round.
    fn round_count(&self, n: usize) -> usize {
        n
    }

    /// Get an array of 1D samples.
    ///
    /// * `n` - The number of samples.
    fn get_1d_array(&mut self, n: usize) -> Vec<Float> {
        self.get_data_mut().get_1d_array(n)
    }

    /// Get an array of 2D samples.
    ///
    /// * `n` - The number of samples.
    fn get_2d_array(&mut self, n: usize) -> Vec<Point2f> {
        self.get_data_mut().get_2d_array(n)
    }

    /// Reset the current sample dimension counter. Returns `true` if `current_pixel_sample_index` <
    /// `samples_per_pixel`; otherwise `false`.
    fn start_next_sample(&mut self) -> bool {
        self.get_data_mut().start_next_sample()
    }

    /// Set the index of the sample in the current pixel to generate next. Returns `true` if
    /// `current_pixel_sample_index` < `samples_per_pixel`; otherwise `false`.
    ///
    /// * `sample_num` - The sample number.
    fn set_sample_number(&mut self, sample_num: usize) -> bool {
        self.get_data_mut().set_sample_number(sample_num)
    }
}

// Implement `Sampler` so `Box<dyn Sampler>` can be passed around where `&dyn Sampler` or `&mut dyn Sampler` can be
// used. This is necessary because we are forced to return heap allocated `Sampler` trait objects from
// `Sampler::clone_sampler()`; at least until a better way is found to clone samplers and pass them around.
impl<S: Sampler + ?Sized> Sampler for Box<S> {
    #[inline]
    fn get_data(&self) -> &SamplerData {
        (**self).get_data()
    }

    #[inline]
    fn get_data_mut(&mut self) -> &mut SamplerData {
        (**self).get_data_mut()
    }

    #[inline]
    fn clone_sampler(&self, seed: u64) -> Box<dyn Sampler> {
        (**self).clone_sampler(seed)
    }

    #[inline]
    fn get_1d(&mut self) -> Float {
        (**self).get_1d()
    }

    #[inline]
    fn get_2d(&mut self) -> Point2f {
        (**self).get_2d()
    }

    #[inline]
    fn start_pixel(&mut self, p: &Point2i) {
        (**self).start_pixel(p);
    }

    #[inline]
    fn get_camera_sample(&mut self, p_raster: &Point2i) -> CameraSample {
        (**self).get_camera_sample(p_raster)
    }

    #[inline]
    fn request_1d_array(&mut self, n: usize) {
        (**self).request_1d_array(n);
    }

    #[inline]
    fn request_2d_array(&mut self, n: usize) {
        (**self).request_2d_array(n);
    }

    #[inline]
    fn round_count(&self, n: usize) -> usize {
        (**self).round_count(n)
    }

    #[inline]
    fn get_1d_array(&mut self, n: usize) -> Vec<Float> {
        (**self).get_1d_array(n)
    }

    #[inline]
    fn get_2d_array(&mut self, n: usize) -> Vec<Point2f> {
        (**self).get_2d_array(n)
    }

    #[inline]
    fn start_next_sample(&mut self) -> bool {
        (**self).start_next_sample()
    }

    #[inline]
    fn set_sample_number(&mut self, sample_num: usize) -> bool {
        (**self).set_sample_number(sample_num)
    }
}

/// Atomic reference counted `Sampler`.
pub type ArcSampler = Arc<dyn Sampler + Send + Sync>;
