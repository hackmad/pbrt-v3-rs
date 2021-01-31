//! Texture Parameters

#![allow(dead_code)]
use super::*;
use crate::core::texture::*;
use std::sync::Arc;

/// Stores texture, geometry and material parameters of different types in hashmaps.
#[derive(Clone)]
pub struct TextureParams {
    /// Floating point textures.
    float_textures: Arc<HashMap<String, ArcTexture<Float>>>,

    /// Spectrum textures.
    spectrum_textures: Arc<HashMap<String, ArcTexture<Spectrum>>>,

    /// Geometry parameters.
    pub geom_params: Arc<ParamSet>,

    /// Material parameters.
    pub mat_params: Arc<ParamSet>,
}

/// Define a macro that can be used to generate a function for finding
/// parameter set item that is stored as a list.
macro_rules! texture_params_find {
    ($func: ident, $t: ty, $paramset_func: ident) => {
        pub fn $func(&mut self, name: &str, mat_default: $t) -> $t {
            let mp = Arc::get_mut(&mut self.mat_params).unwrap();
            let default = mp.$paramset_func(name, mat_default);

            let gp = Arc::get_mut(&mut self.geom_params).unwrap();
            gp.$paramset_func(name, default)
        }
    };
}

impl TextureParams {
    /// Create a new `TextureParams`.
    ///
    /// * `geom_params`        - Geometry parameters.
    /// * `mat_params`        - Material parameters.
    /// * `float_textures`    - Floating point textures.
    /// * `spectrum_textures` - Spectrum textures.
    pub fn new(
        geom_params: Arc<ParamSet>,
        mat_params: Arc<ParamSet>,
        float_textures: Arc<HashMap<String, ArcTexture<Float>>>,
        spectrum_textures: Arc<HashMap<String, ArcTexture<Spectrum>>>,
    ) -> Self {
        Self {
            float_textures: float_textures.clone(),
            spectrum_textures: spectrum_textures.clone(),
            geom_params: geom_params.clone(),
            mat_params: mat_params.clone(),
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
            .map_or(default.clone(), |v| v.clone())
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
            .map_or(default.clone(), |v| v.clone())
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

    /// Prints a report of the material parameter set items that have not been
    /// read after being stored.
    fn report_unused_mat_params<P: fmt::Display>(
        mat: &HashMap<String, ParamSetItem<P>>,
        geom: &HashMap<String, ParamSetItem<P>>,
    ) {
        for (name, param) in mat.iter() {
            if param.looked_up {
                continue;
            }

            // Don't complain about any unused material parameters if their
            // values were provided by a shape parameter.
            if geom.keys().filter(|gp_name| name.eq(*gp_name)).count() > 0 {
                eprintln!("Texture parameter {:} not used.", name);
            }
        }
    }

    /// Prints a report of the geometry and material parameter set items that
    /// have not been read after being stored.
    pub fn report_unused(&self) {
        self.geom_params.report_unused();
        Self::report_unused_mat_params(&self.mat_params.ints, &self.geom_params.ints);
        Self::report_unused_mat_params(&self.mat_params.bools, &self.geom_params.bools);
        Self::report_unused_mat_params(&self.mat_params.floats, &self.geom_params.floats);
        Self::report_unused_mat_params(&self.mat_params.point2fs, &self.geom_params.point2fs);
        Self::report_unused_mat_params(&self.mat_params.vector2fs, &self.geom_params.vector2fs);
        Self::report_unused_mat_params(&self.mat_params.point3fs, &self.geom_params.point3fs);
        Self::report_unused_mat_params(&self.mat_params.vector3fs, &self.geom_params.vector3fs);
        Self::report_unused_mat_params(&self.mat_params.normal3fs, &self.geom_params.normal3fs);
        Self::report_unused_mat_params(&self.mat_params.spectra, &self.geom_params.spectra);
        Self::report_unused_mat_params(&self.mat_params.strings, &self.geom_params.strings);
        Self::report_unused_mat_params(&self.mat_params.textures, &self.geom_params.textures);
    }
}

impl Default for TextureParams {
    /// Initializes a new `TextureParams` with default values.
    fn default() -> Self {
        Self::new(
            Arc::new(ParamSet::new()),
            Arc::new(ParamSet::new()),
            Arc::new(HashMap::new()),
            Arc::new(HashMap::new()),
        )
    }
}
