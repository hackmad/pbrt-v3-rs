//! Parameter Sets

use crate::fileutil::*;
use crate::float_file::parse_float_file;
use crate::geometry::*;
use crate::pbrt::*;
use crate::spectrum::*;
use crate::texture::*;
use std::collections::HashMap;
use std::fmt;

mod paramset_item;
mod texture_params;

// Re-export
pub use paramset_item::*;
pub use texture_params::*;

/// A hashmap of parameter sets stored by name.
pub type ParamSetMap<T> = HashMap<String, ParamSetItem<T>>;

/// Stores parameter set items of different types in hashmaps.
#[derive(Clone)]
pub struct ParamSet {
    pub bools: ParamSetMap<bool>,
    pub ints: ParamSetMap<Int>,
    pub floats: ParamSetMap<Float>,
    pub point2fs: ParamSetMap<Point2f>,
    pub vector2fs: ParamSetMap<Vector2f>,
    pub point3fs: ParamSetMap<Point3f>,
    pub vector3fs: ParamSetMap<Vector3f>,
    pub normal3fs: ParamSetMap<Normal3f>,
    pub spectra: ParamSetMap<Spectrum>,
    pub strings: ParamSetMap<String>,
    pub textures: ParamSetMap<String>,
    pub cached_spectra: HashMap<String, Spectrum>,
}

/// Define a macro that can be used to generate a function for adding/replacing
/// parameter set item.
macro_rules! paramset_add {
    ($func: ident, $t: ty, $paramset: ident) => {
        pub fn $func(&mut self, name: &str, values: &[$t]) {
            let n = String::from(name);
            self.$paramset.insert(n, ParamSetItem::new(values.to_vec()));
        }
    };
}

/// Define a macro that can be used to generate a function for removing
/// parameter set item.
macro_rules! paramset_erase {
    ($func: ident, $paramset: ident) => {
        pub fn $func(&mut self, name: &str) -> bool {
            let n = String::from(name);
            self.$paramset.remove(&n).is_some()
        }
    };
}

/// Define a macro that can be used to generate a function for finding
/// parameter set item that is stored as a single item.
macro_rules! paramset_find_one {
    ($func: ident, $t: ty, $paramset: ident) => {
        pub fn $func(&self, name: &str, default: $t) -> $t {
            let n = String::from(name);
            match self.$paramset.get(&n) {
                Some(param) => {
                    if param.values.len() == 1 {
                        param.values[0].clone()
                    } else {
                        default.clone()
                    }
                }
                None => default.clone(),
            }
        }
    };
}

/// Define a macro that can be used to generate a function for finding
/// parameter set item that is stored as a list.
macro_rules! paramset_find {
    ($func: ident, $t: ty, $paramset: ident) => {
        pub fn $func(&self, name: &str) -> Vec<$t> {
            let n = String::from(name);
            match self.$paramset.get(&n) {
                Some(param) => param.values.clone(),
                None => vec![],
            }
        }
    };
}

/// Define a macro that can be used to print parameter set items.
macro_rules! display_param {
    ($params: expr, $param_type: literal, $formatter: expr) => {
        for (name, param) in $params.iter() {
            let n = param.values.len();
            write!($formatter, "\n\"{} {}\" [", $param_type, name)?;
            if n == 0 {
                write!($formatter, "")?;
            } else {
                write!($formatter, "\n    ")?;

                let mut nc = 0;
                for (i, v) in param.values.iter().enumerate() {
                    let s = if i < n - 1 {
                        format!("{} ", v)
                    } else {
                        format!("{}", v)
                    };

                    nc += s.len();
                    if nc > 80 {
                        write!($formatter, "{}\n    ", s)?;
                        nc = 0;
                    } else {
                        write!($formatter, "{}", s)?;
                    }
                }
                writeln!($formatter, "")?;
            }
            writeln!($formatter, "]")?;
        }
    };
}

impl ParamSet {
    /// Returns a new `ParamSet`.
    pub fn new() -> Self {
        Self {
            bools: HashMap::new(),
            ints: HashMap::new(),
            floats: HashMap::new(),
            point2fs: HashMap::new(),
            vector2fs: HashMap::new(),
            point3fs: HashMap::new(),
            vector3fs: HashMap::new(),
            normal3fs: HashMap::new(),
            spectra: HashMap::new(),
            strings: HashMap::new(),
            textures: HashMap::new(),
            cached_spectra: HashMap::new(),
        }
    }

    paramset_erase!(erase_int, ints);
    paramset_find_one!(find_one_int, Int, ints);
    paramset_find!(find_int, Int, ints);
    paramset_add!(add_int, Int, ints);

    paramset_erase!(erase_bool, bools);
    paramset_find_one!(find_one_bool, bool, bools);
    paramset_find!(find_bool, bool, bools);
    paramset_add!(add_bool, bool, bools);

    paramset_erase!(erase_float, floats);
    paramset_find_one!(find_one_float, Float, floats);
    paramset_find!(find_float, Float, floats);
    paramset_add!(add_float, Float, floats);

    paramset_erase!(erase_point2f, point2fs);
    paramset_find_one!(find_one_point2f, Point2f, point2fs);
    paramset_find!(find_point2f, Point2f, point2fs);
    paramset_add!(add_point2f, Point2f, point2fs);

    paramset_erase!(erase_vector2f, vector2fs);
    paramset_find_one!(find_one_vector2f, Vector2f, vector2fs);
    paramset_find!(find_vector2f, Vector2f, vector2fs);
    paramset_add!(add_vector2f, Vector2f, vector2fs);

    paramset_erase!(erase_point3f, point3fs);
    paramset_find_one!(find_one_point3f, Point3f, point3fs);
    paramset_find!(find_point3f, Point3f, point3fs);
    paramset_add!(add_point3f, Point3f, point3fs);

    paramset_erase!(erase_vector3f, vector3fs);
    paramset_find_one!(find_one_vector3f, Vector3f, vector3fs);
    paramset_find!(find_vector3f, Vector3f, vector3fs);
    paramset_add!(add_vector3f, Vector3f, vector3fs);

    paramset_erase!(erase_normal3f, normal3fs);
    paramset_find_one!(find_one_normal3f, Normal3f, normal3fs);
    paramset_find!(find_normal3f, Normal3f, normal3fs);
    paramset_add!(add_normal3f, Normal3f, normal3fs);

    paramset_erase!(erase_string, strings);
    paramset_find_one!(find_one_string, String, strings);
    paramset_find!(find_string, String, strings);
    paramset_add!(add_string, String, strings);

    paramset_erase!(erase_texture, textures);
    paramset_find_one!(find_one_texture, String, textures);
    paramset_find!(find_texture, String, textures);
    paramset_add!(add_texture, String, textures);

    paramset_erase!(erase_spectrum, spectra);
    paramset_find_one!(find_one_spectrum, Spectrum, spectra);
    paramset_find!(find_spectrum, Spectrum, spectra);

    /// Add/replace an RGB spectrum.
    ///
    /// * `name`   - Parameter name.
    /// * `values` - RGB values in a linear slice.
    pub fn add_rgb_spectrum(&mut self, name: &str, values: &[Float]) {
        let n = values.len();
        assert!(n % 3 == 0, "RGB spectrum values % 3 != 0");

        self.spectra.insert(
            String::from(name),
            ParamSetItem::new(
                (0..n)
                    .step_by(3)
                    .map(|i| Spectrum::from_rgb(&[values[i], values[i + 1], values[i + 2]], None))
                    .collect(),
            ),
        );
    }

    /// Add/replace an XYZ spectrum.
    ///
    /// * `name`   - Parameter name.
    /// * `values` - XYZ values in a linear slice.
    pub fn add_xyz_spectrum(&mut self, name: &str, values: &[Float]) {
        let n = values.len();
        assert!(n % 3 == 0, "XYZ spectrum values % 3 != 0");

        self.spectra.insert(
            String::from(name),
            ParamSetItem::new(
                (0..n)
                    .step_by(3)
                    .map(|i| Spectrum::from_xyz(&[values[i], values[i + 1], values[i + 2]], None))
                    .collect(),
            ),
        );
    }

    /// Add/replace a blackbody spectrum.
    ///
    /// * `name`   - Parameter name.
    /// * `values` - List of (temperature (Kelvin), scale) values in a linear array.
    pub fn add_blackbody_spectrum(&mut self, name: &str, values: &[Float]) {
        let mut n = values.len();
        assert!(n % 2 == 0, "Blackbody spectrum values % 2 != 0");
        n /= 2;

        let lambda = CIE::lambda();
        let mut spectra: Vec<Spectrum> = Vec::with_capacity(n);

        for i in 0..n {
            let le = blackbody_normalized(&lambda, values[2 * i]);
            let samples: Vec<Sample> = lambda
                .iter()
                .zip(le.iter())
                .map(|(l, v)| Sample::new(*l, *v))
                .collect();
            spectra.push(values[2 * i + 1] * Spectrum::from(&samples));
        }

        self.spectra
            .insert(String::from(name), ParamSetItem::new(spectra));
    }

    /// Add/replace a sampled spectrum.
    ///
    /// * `name`   - Parameter name.
    /// * `values` - List of (wavelength, sample) values in a linear array.
    pub fn add_sampled_spectrum(&mut self, name: &str, values: &[Float]) {
        let samples = Sample::list(values);
        let spectra = vec![Spectrum::from(&samples)];
        self.spectra
            .insert(String::from(name), ParamSetItem::new(spectra));
    }

    /// Add/replace a spectra from files.
    ///
    /// * `name`  - Parameter name.
    /// * `paths` - List of paths to the data files.
    /// * `cwd`   - Current working directory for relative path resolution.
    pub fn add_sampled_spectrum_files(&mut self, name: &str, paths: &[String], cwd: &str) {
        let mut spectra: Vec<Spectrum> = vec![];

        for path in paths {
            let p = if is_relative_path(path) && !cwd.is_empty() {
                cwd.to_string() + "/" + &path
            } else {
                path.to_string()
            };

            match absolute_path(&p) {
                Ok(abs_path) => {
                    if let Some(spectrum) = self.cached_spectra.get(&abs_path) {
                        spectra.push(*spectrum);
                        continue;
                    }

                    let spectrum = match parse_float_file(&abs_path) {
                        Ok(vals) => {
                            if vals.len() % 2 > 0 {
                                warn!(
                                    "Extra value found in spectrum file '{}'. Ignoring it.",
                                    path
                                );
                            }
                            let mut samples: Vec<Sample> = Vec::with_capacity(vals.len() / 2);
                            for j in 0..vals.len() / 2 {
                                samples.push(Sample::new(vals[2 * j], vals[2 * j + 1]));
                            }
                            Spectrum::from(&samples)
                        }
                        Err(e) => {
                            warn!(
                                "Unable to read SPD file '{}'. Using black distribution. {}",
                                path, e
                            );
                            Spectrum::ZERO
                        }
                    };
                    self.cached_spectra.insert(abs_path, spectrum);
                    spectra.push(spectrum);
                }
                Err(err) => {
                    error!(
                        "Error reading {}. Using black distribution.\n{}.",
                        path, err
                    );
                    spectra.push(Spectrum::ZERO);
                }
            }
        }

        self.spectra
            .insert(String::from(name), ParamSetItem::new(spectra));
    }

    /// Finds a filename and returns the absolute path to the file.
    ///
    /// * `name` - Parameter name.
    /// * `cwd`  - Current working directory used for relative path handling. If
    ///            it is None, assume path is relative to program path.
    pub fn find_one_filename(&self, name: &str, cwd: Option<&str>) -> Option<String> {
        let mut filename = self.find_one_string(name, String::from(""));
        if filename.is_empty() {
            return None;
        }
        if is_relative_path(&filename) {
            if let Some(c) = cwd {
                if !c.is_empty() {
                    // Path is relative to the parent path of the file being parsed.
                    filename = c.to_string() + "/" + &filename;
                }
            }
        }
        match absolute_path(&filename) {
            Ok(s) => Some(s),
            Err(_) => None,
        }
    }

    /// Clear all parameter set items.
    pub fn clear(&mut self) {
        self.bools.clear();
        self.ints.clear();
        self.floats.clear();
        self.point2fs.clear();
        self.vector2fs.clear();
        self.point3fs.clear();
        self.vector3fs.clear();
        self.normal3fs.clear();
        self.spectra.clear();
        self.strings.clear();
        self.textures.clear();
        self.cached_spectra.clear();
    }
}

impl Default for ParamSet {
    /// Returns the "default value" for `ParamSet`.
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for ParamSet {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        display_param!(self.bools, "bool", f);
        display_param!(self.ints, "integer", f);
        display_param!(self.floats, "float", f);
        display_param!(self.point2fs, "point2", f);
        display_param!(self.vector2fs, "vector2", f);
        display_param!(self.point3fs, "point3", f);
        display_param!(self.vector3fs, "vector3", f);
        display_param!(self.normal3fs, "normal", f);
        display_param!(self.spectra, "color", f);
        display_param!(self.strings, "string", f);
        display_param!(self.textures, "texture", f);
        Ok(())
    }
}
