//! Graphics State

#![allow(dead_code)]
use super::MaterialInstance;
use crate::core::material::*;
use crate::core::medium::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;
use crate::core::texture::*;
use crate::materials::*;
use std::collections::HashMap;
use std::sync::Arc;

lazy_static! {
    /// The global graphics state.
    pub static ref GRAPHICS_STATE: GraphicsState = GraphicsState::default();
}

/// Map of floating point textures.
type FloatTextureMap = HashMap<String, ArcTexture<Float>>;

/// Map of spectrum textures.
type SpectrumTextureMap = HashMap<String, ArcTexture<Spectrum>>;

/// Map of named material instances.
type NamedMaterialMap = HashMap<String, Arc<MaterialInstance>>;

/// Used as a stack to perform hierarchical state management.
#[derive(Clone)]
pub struct GraphicsState {
    /// Tracks the inside medium at a surface boundary.
    pub current_inside_medium: Option<String>,

    /// Tracks the outside medium at a surface boundary.
    pub current_outside_medium: Option<String>,

    /// Stores floating point textures.
    pub float_textures: Arc<FloatTextureMap>,

    /// Indicates if the current floating point texture is shared.
    pub float_textures_shared: bool,

    /// Stores spectrum textures.
    pub spectrum_textures: Arc<SpectrumTextureMap>,

    /// Indicates if the current spectrum texture is shared.
    pub spectrum_textures_shared: bool,

    /// Stores named material instances.
    pub named_materials: Arc<NamedMaterialMap>,

    /// Indicates if the current name material instances is shared.
    pub named_materials_shared: bool,

    /// Current material.
    pub current_material: Option<Arc<MaterialInstance>>,

    /// Current area light parameters.
    pub area_light_params: Arc<ParamSet>,

    /// Current area light name.
    pub area_light: Option<String>,

    /// Reverse surface normal direction for current shape/material.
    pub reverse_orientation: bool,
}

impl GraphicsState {
    /// Returns a material for given shape parameters.
    ///
    /// * `geom_params` - Shape parameters.
    pub fn get_material_for_shape(&self, geom_params: Arc<ParamSet>) -> Option<ArcMaterial> {
        let current_material = self
            .current_material
            .as_ref()
            .expect("GraphicsState has no current material");

        if shape_may_set_material_parameters(geom_params.clone()) {
            // Only create a unique material for the shape if the shape's
            // parameters are (apparently) going to provide values for some of
            // the material parameters.
            let mut mp = TextureParams::new(
                geom_params.clone(),
                current_material.params.clone(),
                self.float_textures.clone(),
                self.spectrum_textures.clone(),
            );
            let mat = make_material(&current_material.name, &mut mp);
            mp.report_unused();
            mat
        } else {
            Some(current_material.material.clone())
        }
    }

    pub fn create_medium_interface(&self) -> MediumInterface {
        todo!()
    }
}

impl Default for GraphicsState {
    /// Initializes a new `GraphicsState`.
    fn default() -> Self {
        // Create a default material.
        let mut mp = TextureParams::default();
        let matte = Arc::new(MatteMaterial::from(&mut mp));
        let current_material = Arc::new(MaterialInstance::new(
            "matte",
            matte,
            Arc::new(ParamSet::new()),
        ));

        Self {
            current_inside_medium: None,
            current_outside_medium: None,
            float_textures: Arc::new(FloatTextureMap::new()),
            float_textures_shared: false,
            spectrum_textures: Arc::new(SpectrumTextureMap::new()),
            spectrum_textures_shared: false,
            named_materials: Arc::new(NamedMaterialMap::new()),
            named_materials_shared: false,
            current_material: Some(current_material),
            area_light_params: Arc::new(ParamSet::new()),
            area_light: None,
            reverse_orientation: false,
        }
    }
}

// Attempt to determine if the ParamSet for a shape may provide a value for
// its material's parameters. Unfortunately, materials don't provide an
// explicit representation of their parameters that we can query and
// cross-reference with the parameter values available from the shape.
//
// Therefore, we'll apply some "heuristics".
fn shape_may_set_material_parameters(ps: Arc<ParamSet>) -> bool {
    for name in ps.textures.keys() {
        // Any texture other than one for an alpha mask is almost certainly
        // for a Material (or is unused!).
        if name != "alpha" && name != "shadowalpha" {
            return true;
        }
    }

    // Special case spheres, which are the most common non-mesh primitive.
    for (name, param) in ps.floats.iter() {
        if param.values.len() == 1 && name != "radius" {
            return true;
        }
    }

    // Extra special case strings, since plymesh uses "filename", curve "type",
    // and loopsubdiv "scheme".
    for (name, param) in ps.strings.iter() {
        if param.values.len() == 1 && name != "filename" && name != "type" && name != "scheme" {
            return true;
        }
    }

    // For all other parameter types, if there is a single value of the
    // parameter, assume it may be for the material. This should be valid
    // (if conservative), since no materials currently take array
    // parameters.
    for param in ps.bools.values() {
        if param.values.len() == 1 {
            return true;
        }
    }

    for param in ps.ints.values() {
        if param.values.len() == 1 {
            return true;
        }
    }

    for param in ps.point2fs.values() {
        if param.values.len() == 1 {
            return true;
        }
    }

    for param in ps.vector2fs.values() {
        if param.values.len() == 1 {
            return true;
        }
    }

    for param in ps.point3fs.values() {
        if param.values.len() == 1 {
            return true;
        }
    }

    for param in ps.vector3fs.values() {
        if param.values.len() == 1 {
            return true;
        }
    }

    for param in ps.normal3fs.values() {
        if param.values.len() == 1 {
            return true;
        }
    }

    for param in ps.spectra.values() {
        if param.values.len() == 1 {
            return true;
        }
    }

    false
}

/// Creates the given type of material from parameter set.
///
/// * `name` - Name.
/// * `mp`   - Parameter set.
fn make_material(name: &str, mp: &mut TextureParams) -> Option<ArcMaterial> {
    match name {
        "matte" => Some(Arc::new(MatteMaterial::from(mp))),
        "plastic" => Some(Arc::new(PlasticMaterial::from(mp))),
        "fourier" => Some(Arc::new(FourierMaterial::from(mp))),
        "mix" => {
            let m1 = mp.find_string("namedmaterial1", String::from(""));
            let m2 = mp.find_string("namedmaterial2", String::from(""));
            let mat1 = match GRAPHICS_STATE.named_materials.get(&m1) {
                None => {
                    eprintln!("Named material '{}' undefined. Using 'matte'.", m1);
                    make_material("matte", mp).unwrap()
                }
                Some(mat) => mat.material.clone(),
            };
            let mat2 = match GRAPHICS_STATE.named_materials.get(&m2) {
                None => {
                    eprintln!("Named material '{}' undefined. Using 'matte'.", m2);
                    make_material("matte", mp).unwrap()
                }
                Some(mat) => mat.material.clone(),
            };

            Some(Arc::new(MixMaterial::from((mp, mat1, mat2))))
        }
        "" => {
            eprintln!("Unable to create material with no name");
            None
        }
        "none" => {
            eprintln!("Unable to create material 'none'.");
            None
        }
        _ => {
            eprintln!("Material '{}' unknown. Using 'matte'.", name);
            Some(Arc::new(MatteMaterial::from(mp)))
        }
    }
}
