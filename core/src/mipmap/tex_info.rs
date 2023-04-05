//! TexInfo

use crate::mipmap::*;
use ordered_float::OrderedFloat;
use std::hash::{Hash, Hasher};

/// Stores the image maps's parameters.
#[derive(Clone)]
pub struct TexInfo {
    /// The path to the image file.
    pub path: String,

    /// Type of filtering to use for mipmaps.
    pub filtering_method: FilteringMethod,

    /// Image wrapping convention.
    pub wrap_mode: ImageWrap,

    /// Scale for the texel values.
    pub scale: Float,

    /// Do gamma correction for the texel values.
    pub gamma: bool,

    /// Used to clamp the ellipse eccentricity (EWA). Set to 0 if EWA is not being used.
    pub max_anisotropy: Float,
}

impl TexInfo {
    /// Create a new `TexInfo`.
    ///
    /// * `path`             - The path to the image file.
    /// * `filtering_method` - Type of filtering to use for mipmaps.
    /// * `wrap_mode`        - Image wrapping convention.
    /// * `scale`            - Scale for the texel values.
    /// * `gamma`            - Do gamma correction for the texel values.
    /// * `max_anisotropy`   - Used to clamp the ellipse eccentricity (EWA). Set to 0 if EWA is not being used.
    pub fn new(
        path: &str,
        filtering_method: FilteringMethod,
        wrap_mode: ImageWrap,
        scale: Float,
        gamma: bool,
        max_anisotropy: Float,
    ) -> Self {
        Self {
            path: String::from(path),
            filtering_method,
            wrap_mode,
            scale,
            gamma,
            max_anisotropy,
        }
    }
}

impl PartialEq for TexInfo {
    /// Checks if all fields are equal.
    ///
    /// * `other` - Another instance of `TexInfo`.
    fn eq(&self, other: &Self) -> bool {
        self.path == other.path
            && self.filtering_method == other.filtering_method
            && self.wrap_mode == other.wrap_mode
            && self.scale == other.scale
            && self.gamma == other.gamma
    }
}

impl Eq for TexInfo {}

impl Hash for TexInfo {
    /// Feeds this value into the given `Hasher`.
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.path.hash(state);
        self.filtering_method.hash(state);
        self.wrap_mode.hash(state);
        self.gamma.hash(state);
        OrderedFloat::from(self.scale).hash(state);
        OrderedFloat::from(self.max_anisotropy).hash(state);
    }
}
