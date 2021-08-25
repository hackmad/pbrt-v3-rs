//! Texture Parameters

use super::*;
use crate::texture::{FloatTextureMap, SpectrumTextureMap};
use std::sync::Arc;

/// Stores texture, geometry and material parameters of different types in hashmaps.
#[derive(Clone)]
pub struct TextureParams {
    /// Floating point textures.
    float_textures: HashMap<String, ArcTexture<Float>>,

    /// Spectrum textures.
    spectrum_textures: HashMap<String, ArcTexture<Spectrum>>,

    /// Geometry parameters.
    pub geom_params: ParamSet,

    /// Material parameters.
    pub mat_params: ParamSet,
}

/// Define a macro that can be used to generate a function for finding
/// parameter set item that is stored as a list.
macro_rules! texture_params_find {
    ($func: ident, $t: ty, $paramset_func: ident) => {
        pub fn $func(&self, name: &str, mat_default: $t) -> $t {
            let default = self.mat_params.$paramset_func(name, mat_default);
            self.geom_params.$paramset_func(name, default)
        }
    };
}

impl TextureParams {
    /// Create a new `TextureParams`.
    ///
    /// * `geom_params`       - Geometry parameters.
    /// * `mat_params`        - Material parameters.
    /// * `float_textures`    - Floating point textures.
    /// * `spectrum_textures` - Spectrum textures.
    pub fn new(
        geom_params: ParamSet,
        mat_params: ParamSet,
        float_textures: FloatTextureMap,
        spectrum_textures: SpectrumTextureMap,
    ) -> Self {
        Self {
            float_textures,
            spectrum_textures,
            geom_params,
            mat_params,
        }
    }

    /// Returns a floating point texture.
    ///
    /// * `name` - Parameter name.
    pub fn get_float_texture(&self, name: &str) -> Option<ArcTexture<Float>> {
        self.float_textures.get(&String::from(name)).cloned()
    }

    /// Returns a floating point texture or a default texture if not found.
    ///
    /// * `name`    - Parameter name.
    /// * `default` - Default texture.
    pub fn get_float_texture_or_else(
        &self,
        name: &str,
        default: ArcTexture<Float>,
    ) -> ArcTexture<Float> {
        self.float_textures
            .get(&String::from(name))
            .map_or(Arc::clone(&default), |v| Arc::clone(&v))
    }

    /// Returns a spectrum point texture.
    ///
    /// * `name` - Parameter name.
    pub fn get_spectrum_texture(&self, name: &str) -> Option<ArcTexture<Spectrum>> {
        self.spectrum_textures.get(&String::from(name)).cloned()
    }

    /// Returns a spectrum point texture or a default texture if not found.
    ///
    /// * `name`    - Parameter name.
    /// * `default` - Default texture.
    pub fn get_spectrum_texture_or_else(
        &self,
        name: &str,
        default: ArcTexture<Spectrum>,
    ) -> ArcTexture<Spectrum> {
        self.spectrum_textures
            .get(&String::from(name))
            .map_or(Arc::clone(&default), |v| Arc::clone(&v))
    }

    texture_params_find!(find_float, Float, find_one_float);
    texture_params_find!(find_string, String, find_one_string);
    texture_params_find!(find_filename, String, find_one_filename);
    texture_params_find!(find_int, Int, find_one_int);
    texture_params_find!(find_bool, bool, find_one_bool);
    texture_params_find!(find_point3f, Point3f, find_one_point3f);
    texture_params_find!(find_vector3f, Vector3f, find_one_vector3f);
    texture_params_find!(find_normal3f, Normal3f, find_one_normal3f);
    texture_params_find!(find_spectrum, Spectrum, find_one_spectrum);
}

impl Default for TextureParams {
    /// Initializes a new `TextureParams` with default values.
    fn default() -> Self {
        Self::new(
            ParamSet::new(),
            ParamSet::new(),
            HashMap::new(),
            HashMap::new(),
        )
    }
}
