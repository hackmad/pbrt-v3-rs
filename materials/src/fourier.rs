//! Fourier Material

use core::bssrdf::*;
use core::interaction::*;
use core::material::*;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::*;
use core::texture::*;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

lazy_static! {
    /// Caches BSDF table data by file path.
    static ref BSDF_TABLES: Mutex<HashMap<String, Arc<FourierBSDFTable>>> = Mutex::new(HashMap::new());
}

/// Implements materials using measured or synthetic BSDF data that has been
/// tabulated into the directional basis.
pub struct FourierMaterial {
    /// Stores the measured Fourier BSDF data.
    bsdf_table: Arc<FourierBSDFTable>,

    /// Bump map.
    bump_map: Option<ArcTexture<Float>>,
}

impl FourierMaterial {
    /// Create a new `FourierMaterial`.
    ///
    /// * `path`     - Path to the Fourier BSDF data file.
    /// * `bump_map` - Optional bump map.
    pub fn new(path: &str, bump_map: Option<ArcTexture<Float>>) -> Self {
        let key = String::from(path);

        // Use preloaded BSDF data if available.
        let mut tables = BSDF_TABLES.lock().unwrap();
        let bsdf_table = if let Some(table) = tables.get(&key) {
            Arc::clone(table)
        } else {
            match FourierBSDFTable::from_file(path) {
                Ok(table) => {
                    let t = Arc::new(table);
                    tables.insert(key, Arc::clone(&t));
                    t
                }
                Err(err) => {
                    panic!("Unable to load file {}. {:}.", path, err);
                }
            }
        };

        Self {
            bsdf_table,
            bump_map,
        }
    }
}

impl Material for FourierMaterial {
    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode.
    /// * `allow_multiple_lobes` - Indicates whether the material should use
    ///                            BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs
    ///                            are available (ignored).
    /// * `bsdf`                 - The computed BSDF.
    /// * `bssrdf`               - The computed BSSSRDF.
    fn compute_scattering_functions<'scene>(
        &self,
        si: &mut SurfaceInteraction<'scene>,
        mode: TransportMode,
        _allow_multiple_lobes: bool,
        bsdf: &mut Option<BSDF>,
        bssrdf: &mut Option<BSSRDF>,
    ) {
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = &self.bump_map {
            Material::bump(self, bump_map, si);
        }

        let mut result = BSDF::new(&si.hit, &si.shading, None);

        // Checking for zero channels works as a proxy for checking whether the
        // table was successfully read from the file.
        if self.bsdf_table.n_channels > 0 {
            result.add(FourierBSDF::new(Arc::clone(&self.bsdf_table), mode));
        }

        *bsdf = Some(result);
        *bssrdf = None;
    }
}

impl From<(&TextureParams, &str)> for FourierMaterial {
    /// Create a Fourier material from given parameter set and current working
    /// directory.
    ///
    /// * `p` - Texture parameter set and current working directory.
    fn from(p: (&TextureParams, &str)) -> Self {
        let (tp, cwd) = p;
        let bump_map = tp.get_float_texture_or_none("bumpmap");
        let path = tp
            .find_filename("bsdffile", Some(cwd))
            .expect("Fourier material bsdffile not specified.");
        Self::new(&path, bump_map)
    }
}
