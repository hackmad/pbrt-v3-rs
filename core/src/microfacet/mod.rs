//! Microfacet Distribution Models

#![allow(dead_code)]
use crate::geometry::*;
use crate::pbrt::*;
use crate::reflection::*;
use std::sync::Arc;

mod beckmann;
mod trowbridge_reitz;

// Re-exports
pub use beckmann::*;
pub use trowbridge_reitz::*;

/// Interface for microfacet distribution models.
pub trait MicrofacetDistribution {
    /// Returns whether or not the visible area is sampled or not.
    fn get_sample_visible_area(&self) -> bool;

    /// Return the differential area of microfacets oriented with the surface
    /// normal `wh`.
    ///
    /// * `wh` - A sample normal from the distrubition of normal vectors.
    fn d(&self, wh: &Vector3f) -> Float;

    /// Returns the invisible masked microfacet area per visible microfacet area.
    ///
    /// * `w` - The direction from camera/viewer.
    fn lambda(&self, w: &Vector3f) -> Float;

    /// Evaluates Smith's masking-shadowing function which gives the fraction of
    /// microfacets that are visible from a given direction.
    ///
    /// * `w` - The direction from camera/viewer.
    fn g1(&self, w: &Vector3f) -> Float {
        1.0 / (1.0 + self.lambda(w))
    }

    /// Returns the fraction of microfacets in a differential area that are
    /// visible from both directions `wo` and `wi`.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    fn g(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        1.0 / (1.0 + self.lambda(wo) + self.lambda(wi))
    }

    /// Returns a sample from the distribution of normal vectors.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    fn sample_wh(&self, wo: &Vector3f, u: &Point2f) -> Vector3f;

    /// Evaluates the PDF for the given outgoing direction and sampled surface
    /// normal.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wh` - A sample normal from the distrubition of normal vectors.
    fn pdf(&self, wo: &Vector3f, wh: &Vector3f) -> Float {
        if self.get_sample_visible_area() {
            self.d(wh) * self.g1(wo) * wo.abs_dot(wh) / abs_cos_theta(wo)
        } else {
            self.d(wh) * abs_cos_theta(wh)
        }
    }
}

/// Atomic reference counted `BSDF`.
pub type ArcMicrofacetDistribution = Arc<dyn MicrofacetDistribution + Send + Sync>;
