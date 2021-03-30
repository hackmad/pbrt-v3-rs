//! Graphics State

#![allow(dead_code)]
use super::material_instance::MaterialInstance;
use super::{
    FloatTextureMap, NamedMaterialMap, SpectrumTextureMap, TransformCache, TransformSet,
    MAX_TRANSFORMS,
};
use crate::accelerators::*;
use crate::cameras::*;
use crate::core::camera::*;
use crate::core::film::*;
use crate::core::filter::*;
use crate::core::geometry::*;
use crate::core::light::*;
use crate::core::material::*;
use crate::core::medium::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;
use crate::core::primitive::*;
use crate::core::sampler::*;
use crate::core::spectrum::*;
use crate::core::texture::*;
use crate::filters::*;
use crate::materials::*;
use crate::samplers::*;
use crate::shapes::*;
use crate::textures::*;
use std::result::Result;
use std::sync::{Arc, Mutex};

/// Used as a stack to perform hierarchical state management.
#[derive(Clone)]
pub struct GraphicsState {
    /// Transform cache.
    pub transform_cache: Arc<Mutex<TransformCache>>,

    /// Tracks the inside medium at a surface boundary.
    pub current_inside_medium: Option<String>,

    /// Tracks the outside medium at a surface boundary.
    pub current_outside_medium: Option<String>,

    /// Stores floating point textures.
    pub float_textures: FloatTextureMap,

    /// Indicates if the current floating point texture is shared.
    pub float_textures_shared: bool,

    /// Stores spectrum textures.
    pub spectrum_textures: SpectrumTextureMap,

    /// Indicates if the current spectrum texture is shared.
    pub spectrum_textures_shared: bool,

    /// Stores named material instances.
    pub named_materials: NamedMaterialMap,

    /// Indicates if the current name material instances is shared.
    pub named_materials_shared: bool,

    /// Current material.
    pub current_material: Option<Arc<MaterialInstance>>,

    /// Current area light parameters.
    pub area_light_params: ParamSet,

    /// Current area light name.
    pub area_light: Option<String>,

    /// Reverse surface normal direction for current shape/material.
    pub reverse_orientation: bool,
}

impl GraphicsState {
    /// Initializes a new `GraphicsState`.
    ///
    /// * `transform_cache` - The `TransformCache`.
    pub fn new(transform_cache: Arc<Mutex<TransformCache>>) -> Self {
        // Create a default material.
        let mp = TextureParams::default();
        let matte = Arc::new(MatteMaterial::from(&mp));
        let current_material = Arc::new(MaterialInstance::new("matte", matte, &ParamSet::new()));

        Self {
            transform_cache: transform_cache.clone(),
            current_inside_medium: None,
            current_outside_medium: None,
            float_textures: FloatTextureMap::new(),
            float_textures_shared: false,
            spectrum_textures: SpectrumTextureMap::new(),
            spectrum_textures_shared: false,
            named_materials: NamedMaterialMap::new(),
            named_materials_shared: false,
            current_material: Some(current_material),
            area_light_params: ParamSet::new(),
            area_light: None,
            reverse_orientation: false,
        }
    }

    /// Returns a material for given shape parameters.
    ///
    /// * `geom_params` - Shape parameters.
    pub fn get_material_for_shape(&self, geom_params: &ParamSet) -> Result<ArcMaterial, String> {
        let current_material = self
            .current_material
            .as_ref()
            .expect("GraphicsState has no current material");

        if self.shape_may_set_material_parameters(geom_params) {
            // Only create a unique material for the shape if the shape's
            // parameters are (apparently) going to provide values for some of
            // the material parameters.
            let mp = TextureParams::new(
                geom_params.clone(),
                current_material.params.clone(),
                self.float_textures.clone(),
                self.spectrum_textures.clone(),
            );
            self.make_material(&current_material.name, &mp)
        } else {
            Ok(current_material.material.clone())
        }
    }

    // Attempt to determine if the ParamSet for a shape may provide a value for
    // its material's parameters. Unfortunately, materials don't provide an
    // explicit representation of their parameters that we can query and
    // cross-reference with the parameter values available from the shape.
    //
    // Therefore, we'll apply some "heuristics".
    fn shape_may_set_material_parameters(&self, ps: &ParamSet) -> bool {
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

    /// Creates the given shape rom parameter set.
    ///
    /// * `name`                - Name.
    /// * `object2world`        - Transformation from object space to world space.
    /// * `world2object`        - Transformation from world space to object space.
    /// * `reverse_orientation` - Indicates whether their surface normal directions.
    /// * `paramset`            - Parameter set.
    pub fn make_shape(
        &self,
        name: &str,
        object2world: ArcTransform,
        world2object: ArcTransform,
        reverse_orientation: bool,
        paramset: &ParamSet,
    ) -> Result<Vec<ArcShape>, String> {
        let p = (paramset, object2world, world2object, reverse_orientation);

        match name {
            "cone" => Ok(vec![Arc::new(Cone::from(p))]),
            "curve" => Ok(Curve::from_props(p)),
            "cylinder" => Ok(vec![Arc::new(Cylinder::from(p))]),
            "disk" => Ok(vec![Arc::new(Disk::from(p))]),
            "hyperboloid" => Ok(vec![Arc::new(Hyperboloid::from(p))]),
            "loopsubdiv" => Ok(LoopSubDiv::from_props(p)),
            "paraboloid" => Ok(vec![Arc::new(Paraboloid::from(p))]),
            "sphere" => Ok(vec![Arc::new(Sphere::from(p))]),
            "trianglemesh" => Ok(TriangleMesh::from_props(p, &self.float_textures)),
            _ => Err(format!("Shape '{}' unknown.", name)),
        }
    }

    /// Creates the given type of material from parameter set.
    ///
    /// * `name` - Name.
    /// * `mp`   - Parameter set.
    pub fn make_material(&self, name: &str, mp: &TextureParams) -> Result<ArcMaterial, String> {
        match name {
            "matte" => Ok(Arc::new(MatteMaterial::from(mp))),
            "plastic" => Ok(Arc::new(PlasticMaterial::from(mp))),
            "fourier" => Ok(Arc::new(FourierMaterial::from(mp))),
            "mix" => {
                let m1 = mp.find_string("namedmaterial1", String::from(""));
                let mat1 = match self.named_materials.get(&m1) {
                    Some(mat) => mat.material.clone(),
                    None => {
                        warn!("Named material '{}' undefined. Using 'matte'.", m1);
                        self.make_material("matte", mp).unwrap()
                    }
                };

                let m2 = mp.find_string("namedmaterial2", String::from(""));
                let mat2 = match self.named_materials.get(&m2) {
                    Some(mat) => mat.material.clone(),
                    None => {
                        warn!("Named material '{}' undefined. Using 'matte'.", m2);
                        self.make_material("matte", mp).unwrap()
                    }
                };

                Ok(Arc::new(MixMaterial::from((mp, mat1, mat2))))
            }
            "" => Err(format!("Unable to create material with no name")),
            "none" => Err(String::from("Unable to create material 'none'.")),
            _ => {
                warn!("Material '{}' unknown. Using 'matte'.", name);
                Ok(Arc::new(MatteMaterial::from(mp)))
            }
        }
    }

    /// Creates a float texture.
    ///
    /// * `name`      - Name.
    /// * `tex2world` - Texture space to world space transform.
    /// * `tp`        - Parameter set.
    pub fn make_float_texture(
        name: &str,
        tex2world: &Transform,
        tp: &TextureParams,
    ) -> Result<ArcTexture<Float>, String> {
        let p = (tp, tex2world);
        match name {
            "bilerp" => Ok(Arc::new(BilerpTexture::<Float>::from(p))),
            "checkerboard" => {
                let dim = p.0.find_int("dimension", 2);
                if dim == 2 {
                    Ok(Arc::new(CheckerboardTexture2D::<Float>::from(p)))
                } else if dim == 3 {
                    Ok(Arc::new(CheckerboardTexture3D::<Float>::from(p)))
                } else {
                    Err(format!(
                        "{} dimensional checkerboard texture not supported",
                        dim
                    ))
                }
            }
            "constant" => Ok(Arc::new(ConstantTexture::<Float>::from(p))),
            "dots" => Ok(Arc::new(DotsTexture::<Float>::from(p))),
            "fbm" => Ok(Arc::new(FBmTexture::<Float>::from(p))),
            "imagemap" => Ok(Arc::new(ImageTexture::<Float>::from(p))),
            "mix" => Ok(Arc::new(MixTexture::<Float>::from(p))),
            "scale" => Ok(Arc::new(ScaleTexture::<Float>::from(p))),
            "windy" => Ok(Arc::new(WindyTexture::<Float>::from(p))),
            _ => Err(format!("Float texture '{}' unknown.", name)),
        }
    }

    /// Creates a spectrum texture.
    ///
    /// * `name`      - Name.
    /// * `tex2world` - Texture space to world space transform.
    /// * `tp`        - Parameter set.
    pub fn make_spectrum_texture(
        name: &str,
        tex2world: &Transform,
        tp: &TextureParams,
    ) -> Result<ArcTexture<Spectrum>, String> {
        let p = (tp, tex2world);
        match name {
            "bilerp" => Ok(Arc::new(BilerpTexture::<Spectrum>::from(p))),
            "checkerboard" => {
                let dim = p.0.find_int("dimension", 2);
                if dim == 2 {
                    Ok(Arc::new(CheckerboardTexture2D::<Spectrum>::from(p)))
                } else if dim == 3 {
                    Ok(Arc::new(CheckerboardTexture3D::<Spectrum>::from(p)))
                } else {
                    Err(format!(
                        "{} dimensional checkerboard texture not supported",
                        dim
                    ))
                }
            }
            "constant" => Ok(Arc::new(ConstantTexture::<Spectrum>::from(p))),
            "dots" => Ok(Arc::new(DotsTexture::<Spectrum>::from(p))),
            "fbm" => Ok(Arc::new(FBmTexture::<Spectrum>::from(p))),
            "imagemap" => Ok(Arc::new(ImageTexture::<Spectrum>::from(p))),
            "marble" => Ok(Arc::new(MarbleTexture::from(p))),
            "mix" => Ok(Arc::new(MixTexture::<Spectrum>::from(p))),
            "scale" => Ok(Arc::new(ScaleTexture::<Spectrum>::from(p))),
            "uv" => Ok(Arc::new(UVTexture::from(p))),
            "windy" => Ok(Arc::new(WindyTexture::<Spectrum>::from(p))),
            _ => Err(format!("Spectrum texture '{}' unknown.", name)),
        }
    }

    /// Creates a medium.
    ///
    /// * `name`         - Name.
    /// * `medium2world` - Medium to world space transform.
    /// * `paramset`     - Parameter set.
    pub fn make_medium(
        _name: &str,
        _medium2world: ArcTransform,
        _paramset: &ParamSet,
    ) -> Result<ArcMedium, String> {
        Err(String::from("GraphicsState::make_medium() not implemented"))
    }

    /// Creates a light.
    ///
    /// * `name`             - Name.
    /// * `light2world`      - Light to world space transform.
    /// * `medium_interface` - Medium interface.
    /// * `paramset`         - Parameter set.
    pub fn make_light(
        _name: &str,
        _light2world: ArcTransform,
        _medium_interface: &MediumInterface,
        _paramset: &ParamSet,
    ) -> Result<ArcLight, String> {
        Err(String::from("GraphicsState::make_light() not implemented"))
    }

    /// Creates an area light.
    ///
    /// NOTE: Upcasting from AreaLight -> Light is not possible. So we return
    /// Result<ArcLight, String>.
    ///
    /// * `name`             - Name.
    /// * `light2world`      - Light to world space transform.
    /// * `medium_interface` - Medium interface.
    /// * `shape`            - Shape
    /// * `paramset`         - Parameter set.
    pub fn make_area_light(
        _name: &str,
        _light2world: ArcTransform,
        _medium_interface: &MediumInterface,
        _shape: ArcShape,
        _paramset: &ParamSet,
    ) -> Result<ArcLight, String> {
        Err(String::from(
            "GraphicsState::make_area_light() not implemented",
        ))
    }

    /// Creates an accelerator.
    ///
    /// * `name`     - Name.
    /// * `prims`    - Primitives.
    /// * `paramset` - Parameter set.
    pub fn make_accelerator(
        name: &str,
        prims: &Vec<ArcPrimitive>,
        paramset: &ParamSet,
    ) -> Result<ArcPrimitive, String> {
        let p = (paramset, prims);
        match name {
            "bvh" => Ok(Arc::new(BVHAccel::from(p))),
            "kdtree" => Ok(Arc::new(KDTreeAccel::from(p))),
            _ => Err(format!("Accelerator '{}' unknown.", name)),
        }
    }

    /// Creates a camera.
    ///
    /// * `name`             - Name.
    /// * `paramset`         - Parameter set.
    /// * `cam2world_set`    - Transform set for camera space to world space
    ///                        transformations.
    /// * `transform_start`  - Start time.
    /// * `transform_end`    - End time.
    /// * `film`             - The film.
    /// * `medium_interface` - The medium interface.
    pub fn make_camera(
        &self,
        name: &str,
        paramset: &ParamSet,
        cam2world_set: &TransformSet,
        transform_start: Float,
        transform_end: Float,
        film: Arc<Film>,
        medium_interface: &MediumInterface,
    ) -> Result<ArcCamera, String> {
        assert!(
            MAX_TRANSFORMS == 2,
            "TransformCache assumes only two transforms"
        );

        let mut transform_cache = self.transform_cache.lock().unwrap();

        let cam2world_start = transform_cache.lookup(cam2world_set[0].clone());
        let cam2world_end = transform_cache.lookup(cam2world_set[1].clone());

        let animated_cam2world = AnimatedTransform::new(
            cam2world_start,
            cam2world_end,
            transform_start,
            transform_end,
        );

        let p = (
            paramset,
            &animated_cam2world,
            film,
            medium_interface.outside.clone(),
        );

        match name {
            "environment" => Ok(Arc::new(EnvironmentCamera::from(p))),
            "orthographic" => Ok(Arc::new(OrthographicCamera::from(p))),
            "perspective" => Ok(Arc::new(PerspectiveCamera::from(p))),
            "realistic" => Ok(Arc::new(RealisticCamera::from(p))),
            _ => Err(format!("Camera '{}' unknown.", name)),
        }
    }

    /// Creates a sampler.
    ///
    /// * `name`     - Name.
    /// * `paramset` - Parameter set.
    /// * `film`     - The film.
    pub fn make_sampler(
        name: &str,
        paramset: &ParamSet,
        film: Arc<Film>,
    ) -> Result<ArcSampler, String> {
        let p = (paramset, film.clone().get_sample_bounds());

        match name {
            "02sequence" => Ok(Arc::new(ZeroTwoSequenceSampler::from(p))),
            "lowdiscrepency" => Ok(Arc::new(ZeroTwoSequenceSampler::from(p))),
            "halton" => Ok(Arc::new(HaltonSampler::from(p))),
            "maxmindist" => Ok(Arc::new(HaltonSampler::from(p))),
            "random" => Ok(Arc::new(RandomSampler::from(p))),
            "sobol" => Ok(Arc::new(SobolSampler::from(p))),
            "stratified" => Ok(Arc::new(StratifiedSampler::from(p))),
            _ => Err(format!("Sampler '{}' unknown.", name)),
        }
    }

    /// Creates a filter.
    ///
    /// * `name`     - Name.
    /// * `paramset` - Parameter set.
    pub fn make_filter(name: &str, paramset: &ParamSet) -> Result<ArcFilter, String> {
        match name {
            "box" => Ok(Arc::new(BoxFilter::from(paramset))),
            "gaussian" => Ok(Arc::new(GaussianFilter::from(paramset))),
            "mitchell" => Ok(Arc::new(MitchellFilter::from(paramset))),
            "sinc" => Ok(Arc::new(LanczosSincFilter::from(paramset))),
            "triangle" => Ok(Arc::new(TriangleFilter::from(paramset))),
            _ => Err(format!("Sampler '{}' unknown.", name)),
        }
    }

    /// Creates a film.
    ///
    /// * `name`     - Name.
    /// * `paramset` - Parameter set.
    /// * `filter`   - Filter.
    pub fn make_film(
        name: &str,
        paramset: &ParamSet,
        filter: ArcFilter,
    ) -> Result<Arc<Film>, String> {
        match name {
            "image" => Ok(Arc::new(Film::from((paramset, filter)))),
            _ => Err(format!("Film '{}' unknown.", name)),
        }
    }
}
