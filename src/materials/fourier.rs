//! Fourier Material

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::material::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;
use crate::core::reflection::*;
use crate::core::texture::*;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

lazy_static! {
    /// Caches BSDF table data by file path.
    static ref BSDF_TABLES: Mutex<HashMap<String, Arc<FourierBSDFTable>>> = Mutex::new(HashMap::new());
}

/// Implements materials using measured or synthetic BSDF data.
pub struct FourierMaterial {
    /// Stores the measured Fourier BSDF data.
    bsdf_table: Arc<FourierBSDFTable>,

    /// Bump map.
    bump_map: Option<ArcTexture<Float>>,
}

impl FourierMaterial {
    /// Create a new `FourierMaterial`.
    ///
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
            bump_map: bump_map.clone(),
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
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        // Perform bump mapping with `bump_map`, if present.
        if let Some(bump_map) = self.bump_map.clone() {
            Material::bump(self, bump_map, si);
        }

        let mut bsdf = BSDF::new(&si.clone(), None);

        // Checking for zero channels works as a proxy for checking whether the
        // table was successfully read from the file.
        if self.bsdf_table.n_channels > 0 {
            bsdf.add(Arc::new(FourierBSDF::new(
                Arc::clone(&self.bsdf_table),
                mode,
            )));
        }

        si.bsdf = Some(Arc::new(bsdf));
    }
}

impl From<&TextureParams> for FourierMaterial {
    /// Create a Fourier material from given parameter set.
    ///
    /// * `tp` - Texture parameter set.
    fn from(tp: &TextureParams) -> Self {
        let bump_map = tp.get_float_texture("bumpmap");
        let path = tp.find_filename("bsdfffile", String::from(""));
        Self::new(&path, bump_map)
    }
}
