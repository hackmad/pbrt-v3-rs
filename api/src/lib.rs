//! The API

#[macro_use]
extern crate log;
#[macro_use]
extern crate pest_derive;

mod graphics_state;
mod material_instance;
mod render_options;
mod transform_cache;
mod transform_set;

use accelerators::*;
use core::geometry::*;
use core::light::*;
use core::medium::*;
use core::paramset::*;
use core::pbrt::*;
use core::primitive::*;
use core::primitives::*;
use graphics_state::*;
use material_instance::*;
use render_options::*;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use transform_cache::*;
use transform_set::*;

pub mod parser;

/// Map of named material instances.
pub type NamedMaterialMap = HashMap<String, Arc<MaterialInstance>>;

/// Enumerations for API state as we parse the PBRT file format.
#[derive(Copy, Clone, PartialEq)]
pub enum ApiState {
    /// Not initialized.
    Uninitialized,
    /// In the Options block.
    OptionsBlock,
    /// In the World block.
    WorldBlock,
}

/// State of the PBRT API.
pub struct Api {
    /// Current state of parsing PBRT file.
    current_api_state: ApiState,

    /// Current set of transformations.
    current_transforms: TransformSet,

    /// Active transformation (start time, end time or both) to which
    /// new transformations will be applied.
    active_transform_bits: usize,

    /// Stores transformations for named coordinate systems.
    named_coordinate_systems: HashMap<String, TransformSet>,

    /// Render options.
    render_options: RenderOptions,

    /// Graphics state.
    graphics_state: GraphicsState,

    /// Used as a stack for the graphics state.
    pushed_graphics_states: Vec<GraphicsState>,

    /// Used as a stack for the transformations.
    pushed_transforms: Vec<TransformSet>,

    /// Used as a stack for the active transformations.
    pushed_active_transform_bits: Vec<usize>,

    /// Caches the transforms.
    transform_cache: Arc<Mutex<TransformCache>>,

    /// Current working directory. Used to resolve relative paths.
    cwd: String,
}

impl Api {
    /// Returns a newly initialized API.
    pub fn new() -> Self {
        let transform_cache = Arc::new(Mutex::new(TransformCache::new()));
        let cwd = std::env::current_dir()
            .unwrap()
            .as_path()
            .to_str()
            .unwrap()
            .to_string();

        Self {
            current_api_state: ApiState::Uninitialized,
            current_transforms: TransformSet::default(),
            active_transform_bits: ALL_TRANSFORM_BITS,
            named_coordinate_systems: HashMap::new(),
            render_options: RenderOptions::new(),
            graphics_state: GraphicsState::new(Arc::clone(&transform_cache), &cwd),
            pushed_graphics_states: vec![],
            pushed_transforms: vec![],
            pushed_active_transform_bits: vec![],
            transform_cache: Arc::clone(&transform_cache),
            cwd,
        }
    }

    /* API Methods */

    /// Set current working directory.
    ///
    /// * `path` - The path.
    pub fn set_current_working_dir(&mut self, path: &str) {
        self.cwd = path.to_string();
        self.graphics_state.set_current_working_dir(path);
    }

    /// API Initialization.
    pub fn pbrt_init(&mut self) {
        if self.current_api_state != ApiState::Uninitialized {
            error!("pbrt_init() has already been called.");
        }
        self.current_api_state = ApiState::OptionsBlock;
    }

    /// API Cleanup.
    pub fn pbrt_cleanup(&mut self) {
        if self.current_api_state == ApiState::Uninitialized {
            error!("pbrt_cleanup() called without pbrt_init().");
        } else if self.current_api_state == ApiState::WorldBlock {
            error!("pbrt_cleanup() called while inside world block.");
        }
        self.current_api_state = ApiState::Uninitialized;
    }

    /// Set current tranformation matrix to the identity matrix.
    pub fn pbrt_identity(&mut self) {
        if self.verify_initialized("Identity") {
            for i in 0..MAX_TRANSFORMS {
                if self.active_transform_bits & (1 << i) > 0 {
                    self.current_transforms[i] = Arc::new(Transform::IDENTITY);
                }
            }
        }
    }

    /// Apply a translation to the active transformation.
    ///
    /// * `dx` - Translation in x-direction.
    /// * `dy` - Translation in y-direction.
    /// * `dz` - Translation in z-direction.
    pub fn pbrt_translate(&mut self, dx: Float, dy: Float, dz: Float) {
        if self.verify_initialized("Translate") {
            let transform = Transform::translate(&Vector3f::new(dx, dy, dz));
            for i in 0..MAX_TRANSFORMS {
                if self.active_transform_bits & (1 << i) > 0 {
                    let t = self.current_transforms[i].as_ref() * &transform;
                    self.current_transforms[i] = Arc::new(t);
                }
            }
        }
    }

    /// Setthe active transformation to a given transformation matrix.
    ///
    /// * `tr` - The transformation matrix in column-major form.
    ///          |0, 4,  8, 12|
    ///          |1, 5,  9, 13|
    ///          |2, 6, 10, 14|
    ///          |3, 7, 11, 15|
    pub fn pbrt_transform(&mut self, tr: &[Float; 16]) {
        if self.verify_initialized("Transform") {
            for i in 0..MAX_TRANSFORMS {
                let t = Transform::from(Matrix4x4::new(
                    tr[0], tr[4], tr[8], tr[12], tr[1], tr[5], tr[9], tr[13], tr[2], tr[6], tr[10],
                    tr[14], tr[3], tr[7], tr[11], tr[15],
                ));
                self.current_transforms[i] = Arc::new(t);
            }
        }
    }

    /// Concatenate a transformation to the active transformation.
    ///
    /// * `tr` - Transformation matrix in column-major format.
    ///          |0, 4,  8, 12|
    ///          |1, 5,  9, 13|
    ///          |2, 6, 10, 14|
    ///          |3, 7, 11, 15|
    pub fn pbrt_concat_transform(&mut self, tr: &[Float; 16]) {
        if self.verify_initialized("ConcatTransform") {
            let transform = Transform::from(Matrix4x4::new(
                tr[0], tr[4], tr[8], tr[12], tr[1], tr[5], tr[9], tr[13], tr[2], tr[6], tr[10],
                tr[14], tr[3], tr[7], tr[11], tr[15],
            ));
            for i in 0..MAX_TRANSFORMS {
                let t = self.current_transforms[i].as_ref() * &transform;
                self.current_transforms[i] = Arc::new(t);
            }
        }
    }

    /// Apply a rotation about a vector to the active transformation.
    ///
    /// * `angle` - Angle of rotation in degrees.
    /// * `dx`    - x-component of axis vector.
    /// * `dy`    - y-component of axis vector.
    /// * `dz`    - z-component of axis vector.
    pub fn pbrt_rotate(&mut self, angle: Float, dx: Float, dy: Float, dz: Float) {
        if self.verify_initialized("Rotate") {
            let transform = Transform::rotate_axis(angle, &Vector3f::new(dx, dy, dz));
            for i in 0..MAX_TRANSFORMS {
                if self.active_transform_bits & (1 << i) > 0 {
                    let t = self.current_transforms[i].as_ref() * &transform;
                    self.current_transforms[i] = Arc::new(t);
                }
            }
        }
    }

    /// Apply a scale to the active transformation.
    ///
    /// * `sx` - Scale factor in x-direction.
    /// * `sy` - Scale factor in y-direction.
    /// * `sz` - Scale factor in z-direction.
    pub fn pbrt_scale(&mut self, sx: Float, sy: Float, sz: Float) {
        if self.verify_initialized("Scale") {
            let transform = Transform::scale(sx, sy, sz);
            for i in 0..MAX_TRANSFORMS {
                if self.active_transform_bits & (1 << i) > 0 {
                    let t = self.current_transforms[i].as_ref() * &transform;
                    self.current_transforms[i] = Arc::new(t);
                }
            }
        }
    }

    /// Apply a transformation to the active transformation, to point the camera
    /// in a given direction.
    ///
    /// * `ex` - Camera x-position.
    /// * `ey` - Camera y-position.
    /// * `ez` - Camera z-position.
    /// * `lx` - Look at point x-position.
    /// * `ly` - Look at point y-position.
    /// * `lz` - Look at point z-position.
    /// * `ux` - Up vector x-component.
    /// * `uy` - Up vector y-component.
    /// * `uz` - Up vector z-component.
    pub fn pbrt_look_at(
        &mut self,
        ex: Float,
        ey: Float,
        ez: Float,
        lx: Float,
        ly: Float,
        lz: Float,
        ux: Float,
        uy: Float,
        uz: Float,
    ) {
        if self.verify_initialized("LookAt") {
            let transform = Transform::look_at(
                &Point3f::new(ex, ey, ez),
                &Point3f::new(lx, ly, lz),
                &Vector3f::new(ux, uy, uz),
            );
            for i in 0..MAX_TRANSFORMS {
                if self.active_transform_bits & (1 << i) > 0 {
                    let t = self.current_transforms[i].as_ref() * &transform;
                    self.current_transforms[i] = Arc::new(t);
                }
            }
        }
    }

    /// Stores the current transformation matrices as named coordinate system.
    ///
    /// * `name` - The coordinate system name.
    pub fn pbrt_coordinate_system(&mut self, name: String) {
        if self.verify_initialized("CoordinateSystem") {
            let transforms = self.current_transforms.clone();
            self.named_coordinate_systems.insert(name, transforms);
        }
    }

    /// Restores the current transformation matrices from a named coordinate system.
    ///
    /// * `name` - The coordinate system name.
    pub fn pbrt_coord_sys_transform(&mut self, name: String) {
        if self.verify_initialized("CoordSysTransform") {
            if let Some(transforms) = self.named_coordinate_systems.get(&name) {
                self.current_transforms = (*transforms).clone();
            } else {
                warn!("Couldn't find named coordinate system '{}'.", name);
            }
        }
    }

    /// Set the active transformations to affect both starting and ending time.
    pub fn pbrt_active_transform_all(&mut self) {
        self.active_transform_bits = ALL_TRANSFORM_BITS;
    }

    /// Set the active transformations to affect ending time only.
    pub fn pbrt_active_transform_end_time(&mut self) {
        self.active_transform_bits = END_TRANSFORM_BITS;
    }

    /// Set the active transformations to affect starting time only.
    pub fn pbrt_active_transform_start_time(&mut self) {
        self.active_transform_bits = START_TRANSFORM_BITS;
    }

    /// Set the time index for the starting and ending of transformations.
    ///
    /// * `start` - Starting time.
    /// * `end`   - Ending time.
    pub fn pbrt_transform_times(&mut self, start: Float, end: Float) {
        if self.verify_options("TransformTimes") {
            self.render_options.transform_start_time = start;
            self.render_options.transform_end_time = end;
        }
    }

    /// Set the filter type and parameters used for the film.
    ///
    /// * `name`   - Filter type name.
    /// * `params` - Filter parameters.
    pub fn pbrt_pixel_filter(&mut self, name: String, params: &ParamSet) {
        if self.verify_options("PixelFilter") {
            self.render_options.filter_name = name;
            self.render_options.filter_params = params.clone();
        }
    }

    /// Set the film type and parameters.
    ///
    /// * `name`   - Film type name.
    /// * `params` - Film parameters.
    pub fn pbrt_film(&mut self, film_type: String, params: &ParamSet) {
        if self.verify_options("Film") {
            self.render_options.film_name = film_type;
            self.render_options.film_params = params.clone();
        }
    }

    /// Set the sampler type and parameters.
    ///
    /// * `name`   - Sampler type name.
    /// * `params` - Sampler parameters.
    pub fn pbrt_sampler(&mut self, name: String, params: &ParamSet) {
        if self.verify_options("Sampler") {
            self.render_options.sampler_name = name;
            self.render_options.sampler_params = params.clone();
        }
    }

    /// Set the accelerator type and parameters.
    ///
    /// * `name`   - Accelerator type name.
    /// * `params` - Accelerator parameters.
    pub fn pbrt_accelerator(&mut self, name: String, params: &ParamSet) {
        if self.verify_options("Accelerator") {
            self.render_options.accelerator_name = name;
            self.render_options.accelerator_params = params.clone();
        }
    }

    /// Set the integrator type and parameters.
    ///
    /// * `name`   - Integrator type name.
    /// * `params` - Integrator parameters.
    pub fn pbrt_integrator(&mut self, name: String, params: &ParamSet) {
        if self.verify_options("Integrator") {
            self.render_options.integrator_name = name;
            self.render_options.integrator_params = params.clone();
        }
    }

    /// Set the camera type and parameters. Also sets the camera-to-world transformation
    /// using the inverse of the current transformation matrices.
    ///
    /// * `name`   - Camera type name.
    /// * `params` - Camera parameters.
    pub fn pbrt_camera(&mut self, name: String, params: &ParamSet) {
        if self.verify_options("Camera") {
            self.render_options.camera_name = name;
            self.render_options.camera_params = params.clone();
            self.render_options.camera_to_world = self.current_transforms.inverse();
        }
    }

    /// Create a named medium for participating media.
    ///
    /// * `name`   - Medium name.
    /// * `params` - Medium parameters.
    pub fn pbrt_make_named_medium(&mut self, name: String, params: &ParamSet) {
        if self.verify_initialized("MakeNamedMedium") {
            self.warn_if_animated_transform("MakeNamedMedium");

            let medium_type = params.find_one_string("type", String::new());
            if medium_type.is_empty() {
                error!("No parameter string 'type' found in MakeNamedMedium.");
            } else if let Ok(medium) =
                GraphicsState::make_medium(&name, self.current_transforms[0].clone(), params)
            {
                self.render_options.named_media.insert(name, medium);
            }
        }
    }

    /// Store the names of 2 participating media that form a boundary at some surface.
    ///
    /// * `inside_name`  - Inside medium name
    /// * `outside_name` - Outside medium name.
    pub fn pbrt_medium_interface(&mut self, inside_name: String, outside_name: String) {
        if self.verify_initialized("MediumInterface") {
            self.graphics_state.current_inside_medium = Some(inside_name);
            self.graphics_state.current_outside_medium = Some(outside_name);
            self.render_options.have_scattering_media = true;
        }
    }

    /// Begin world description.
    pub fn pbrt_world_begin(&mut self) {
        if self.verify_options("WorldBegin") {
            self.current_api_state = ApiState::WorldBlock;
            for i in 0..MAX_TRANSFORMS {
                self.current_transforms[i] = Arc::new(Transform::IDENTITY);
            }
            self.active_transform_bits = ALL_TRANSFORM_BITS;
            self.named_coordinate_systems
                .insert(String::from("world"), self.current_transforms.clone());
        }
    }

    /// End world description.
    pub fn pbrt_world_end(&mut self) {
        if self.verify_world("WorldEnd") {
            // Ensure there are no pushed graphics states.
            while !self.pushed_graphics_states.is_empty() {
                warn!("Missing end to pbrtAttributeBegin().");
                self.pushed_graphics_states.pop();
                self.pushed_transforms.pop();
            }

            while !self.pushed_transforms.is_empty() {
                warn!("Missing end to pbrtTransformBegin().");
                self.pushed_transforms.pop();
            }

            // Create scene and render.
            let integrator = match self.render_options.make_integrator(&self.graphics_state) {
                Ok(integrator) => integrator,
                Err(err) => panic!("Error creating integrator. {}", err),
            };

            let scene = self.render_options.make_scene();
            integrator.render(&scene);

            // Clean up after rendering.
            let mut transform_cache = self.transform_cache.lock().unwrap();
            transform_cache.clear();

            self.graphics_state = GraphicsState::new(Arc::clone(&self.transform_cache), &self.cwd);
            self.current_api_state = ApiState::OptionsBlock;
            self.current_transforms.reset();

            self.active_transform_bits = ALL_TRANSFORM_BITS;
            self.named_coordinate_systems.clear();

            // TODO Clear image texture caches for float and spectrum textures
            // once we add this functionality to crate::textures::image_map
        }
    }

    /// Begin an attribute section where current graphics state can be
    /// pushed onto the stack.
    pub fn pbrt_attribute_begin(&mut self) {
        if self.verify_world("AttributeBegin") {
            self.pushed_graphics_states
                .push(self.graphics_state.clone());
            self.graphics_state.float_textures_shared = true;
            self.graphics_state.spectrum_textures_shared = true;
            self.graphics_state.named_materials_shared = true;
            self.pushed_transforms.push(self.current_transforms.clone());
            self.pushed_active_transform_bits
                .push(self.active_transform_bits);
        }
    }

    /// End the attribute section where current graphics state can be
    /// popped off the stack and restored.
    pub fn pbrt_attribute_end(&mut self) {
        if self.verify_world("AttributeEnd") {
            if let Some(graphics_state) = self.pushed_graphics_states.pop() {
                self.graphics_state = graphics_state;
            } else {
                error!("Unmatched pbrtAttributeEnd() encountered. Ignoring it.");
            }
            if let Some(transforms) = self.pushed_transforms.pop() {
                self.current_transforms = transforms;
            }
            if let Some(active_transform_bits) = self.pushed_active_transform_bits.pop() {
                self.active_transform_bits = active_transform_bits;
            }
        }
    }

    /// Save the transformation matrix on the stack independantly of the
    /// graphics state.
    pub fn pbrt_transform_begin(&mut self) {
        if self.verify_world("TransformBegin") {
            self.pushed_transforms.push(self.current_transforms.clone());
            self.pushed_active_transform_bits
                .push(self.active_transform_bits);
        }
    }

    /// Restore the transformation matrix off the stack.
    pub fn pbrt_transform_end(&mut self) {
        if self.verify_world("TransformEnd") {
            if let Some(transforms) = self.pushed_transforms.pop() {
                self.current_transforms = transforms;
            } else {
                error!("Unmatched pbrtTransformEnd() encountered. Ignoring it.");
            }
            if let Some(active_transform_bits) = self.pushed_active_transform_bits.pop() {
                self.active_transform_bits = active_transform_bits;
            }
        }
    }

    /// Create a named texture of particular type.
    ///
    /// * `name`         - Texture name for storage.
    /// * `texture_type` - Texture type (float or color or spectrum).
    /// * `tex_name`     - Texture name (bilerp, checkerboard, etc).
    /// * `params`       - Texture parameters.
    pub fn pbrt_texture(
        &mut self,
        name: String,
        tex_type: String,
        tex_class: String,
        params: &ParamSet,
    ) {
        if self.verify_world("Texture") {
            let tp = TextureParams::new(
                params.clone(),
                params.clone(),
                self.graphics_state.float_textures.clone(),
                self.graphics_state.spectrum_textures.clone(),
            );

            if tex_type == "float" {
                // Create `Float` texture and store in `float_textures`.
                if self.graphics_state.float_textures.contains_key(&name) {
                    warn!("Texture '{}' being redefined.", name);
                }

                self.warn_if_animated_transform("Texture");

                if let Ok(ft) = self.graphics_state.make_float_texture(
                    &tex_class,
                    Arc::clone(&self.current_transforms[0]),
                    &tp,
                ) {
                    if self.graphics_state.float_textures_shared {
                        let ftm = self.graphics_state.float_textures.clone();
                        self.graphics_state.float_textures = ftm;
                        self.graphics_state.float_textures_shared = false;
                    }
                    self.graphics_state.float_textures.insert(name, ft);
                }
            } else if tex_type == "color" || tex_type == "spectrum" {
                // Create `colour` texture and store in `spectrum_textures`.
                if self.graphics_state.spectrum_textures.contains_key(&name) {
                    warn!("Texture '{}' being redefined.", name);
                }

                self.warn_if_animated_transform("Texture");

                if let Ok(st) = self.graphics_state.make_spectrum_texture(
                    &tex_class,
                    Arc::clone(&self.current_transforms[0]),
                    &tp,
                ) {
                    if self.graphics_state.spectrum_textures_shared {
                        let stm = self.graphics_state.spectrum_textures.clone();
                        self.graphics_state.spectrum_textures = stm;
                        self.graphics_state.spectrum_textures_shared = false;
                    }
                    self.graphics_state.spectrum_textures.insert(name, st);
                }
            } else {
                error!("Texture type '{}' unknown.", tex_type);
            }
        }
    }

    /// Specify the current material type and parameters.
    ///
    /// * `name`   - Material type (matte, fourier, etc).
    /// * `params` - Material parameters.
    pub fn pbrt_material(&mut self, name: String, params: &ParamSet) {
        if self.verify_world("Material") {
            let empty_params = ParamSet::new();
            let mp = TextureParams::new(
                params.clone(),
                empty_params,
                self.graphics_state.float_textures.clone(),
                self.graphics_state.spectrum_textures.clone(),
            );
            if let Ok(mtl) = self.graphics_state.make_material(&name, &mp) {
                self.graphics_state.current_material = Some(Arc::new(MaterialInstance::new(
                    &name,
                    Arc::clone(&mtl),
                    params,
                )))
            }
        }
    }

    /// Create a named material with the given parameters.
    ///
    /// * `name`   - Material name.
    /// * `params` - Material parameters.
    pub fn pbrt_make_named_material(&mut self, name: String, params: &ParamSet) {
        if self.verify_world("MakeNamedMaterial") {
            let empty_params = ParamSet::new();
            let mp = TextureParams::new(
                params.clone(),
                empty_params,
                self.graphics_state.float_textures.clone(),
                self.graphics_state.spectrum_textures.clone(),
            );

            self.warn_if_animated_transform("MakeNamedMaterial");

            let mat_name = mp.find_string("type", String::new());
            if mat_name.is_empty() {
                error!("No parameter string 'type' found in MakeNamedMaterial.");
            } else if let Ok(mtl) = self.graphics_state.make_material(&mat_name, &mp) {
                if self.graphics_state.named_materials.contains_key(&name) {
                    warn!("Named material '{}' redefined.", name);
                }
                if self.graphics_state.named_materials_shared {
                    let nm = self.graphics_state.named_materials.clone();
                    self.graphics_state.named_materials = nm;
                    self.graphics_state.named_materials_shared = false;
                }
                let mtli = Arc::new(MaterialInstance::new(&mat_name, Arc::clone(&mtl), params));
                self.graphics_state.named_materials.insert(name, mtli);
            }
        }
    }

    /// Set a named material as current material.
    ///
    /// * `name`   - Material name.
    pub fn pbrt_named_material(&mut self, name: String) {
        if self.verify_world("NamedMaterial") {
            if let Some(mtl) = self.graphics_state.named_materials.get(&name) {
                self.graphics_state.current_material = Some((*mtl).clone());
            } else {
                error!("NamedMaterial '{}' unknown.", name);
            }
        }
    }

    /// Define a light source.
    ///
    /// * `name`   - Light type (point, spot, etc)
    /// * `params` - Light parameters.
    pub fn pbrt_light_source(&mut self, name: String, params: &ParamSet) {
        if self.verify_world("LightSource") {
            self.warn_if_animated_transform("LightSource");

            let mi = self.create_medium_interface();
            let light2world = self.current_transforms[0].clone();
            match self
                .graphics_state
                .make_light(&name, light2world, &mi, params)
            {
                Ok(lt) => self.render_options.lights.push(lt),
                Err(err) => error!("{}", err),
            }
        }
    }

    /// Define an area light source.
    ///
    /// * `name`   - Area light name.
    /// * `params` - Area light parameters.
    pub fn pbrt_area_light_source(&mut self, name: String, params: &ParamSet) {
        if self.verify_world("AreaLightSource") {
            self.graphics_state.area_light = Some(name);
            self.graphics_state.area_light_params = params.clone();
        }
    }

    /// Define a shape.
    ///
    /// * `name`   - Shape type (e.g. sphere, cone, etc)
    /// * `params` - Shape parameters.
    pub fn pbrt_shape(&mut self, name: String, params: &ParamSet) {
        if self.verify_world("Shape") {
            let mut prims: Vec<ArcPrimitive> = vec![];
            let mut area_lights: Vec<ArcLight> = vec![]; // Upcasting AreaLight -> Light not possible.

            if !self.current_transforms.is_animated() {
                // Initialize `prims` and `area_lights` for static shape.

                // Create shapes for shape `name`.
                let mut transform_cache = self.transform_cache.lock().unwrap();
                let tr = self.current_transforms[0].clone();
                let tr_inv = tr.inverse();
                let obj2world = transform_cache.lookup(&tr);
                let world2obj = transform_cache.lookup(&tr_inv);
                let shapes = self
                    .graphics_state
                    .make_shape(
                        &name,
                        Arc::clone(&obj2world),
                        Arc::clone(&world2obj),
                        self.graphics_state.reverse_orientation,
                        params,
                    )
                    .unwrap();

                if shapes.is_empty() {
                    return;
                }

                let mtl = self.graphics_state.get_material_for_shape(params).unwrap();
                let mi = self.create_medium_interface();

                for shape in shapes.iter() {
                    // Possibly create area light for shape.
                    if let Some(area_light) = self.graphics_state.area_light.clone() {
                        if let Ok(area) = GraphicsState::make_area_light(
                            &area_light,
                            self.current_transforms[0].clone(),
                            &mi,
                            Arc::clone(shape),
                            params,
                        ) {
                            area_lights.push(area);
                        }
                    }

                    let prim = GeometricPrimitive::new(
                        Arc::clone(shape),
                        Arc::clone(&mtl),
                        None,
                        mi.clone(),
                    );
                    prims.push(Arc::new(prim));
                }
            } else {
                // Initialize `prims` and `area_lights` for animated shape.

                // Create initial shape or shapes for animated shape.
                if self.graphics_state.area_light.is_none() {
                    warn!("Ignoring currently set area light when creating 'animated shape'.");
                }

                let mut transform_cache = self.transform_cache.lock().unwrap();
                let identity = transform_cache.lookup(&Transform::IDENTITY);
                let shapes = self
                    .graphics_state
                    .make_shape(
                        &name,
                        Arc::clone(&identity),
                        Arc::clone(&identity),
                        self.graphics_state.reverse_orientation,
                        params,
                    )
                    .unwrap();

                if shapes.is_empty() {
                    return;
                }

                // Create `GeometricPrimitive`(s) for animated shape.
                let mtl = self.graphics_state.get_material_for_shape(params).unwrap();
                let mi = self.create_medium_interface();

                for shape in shapes.iter() {
                    let prim = GeometricPrimitive::new(
                        Arc::clone(shape),
                        Arc::clone(&mtl),
                        None,
                        mi.clone(),
                    );
                    prims.push(Arc::new(prim));
                }

                // Create single `TransformedPrimitive` for `prims`.

                // Get `animated_object_to_world` transform for shape.
                let obj2world = [
                    transform_cache.lookup(&self.current_transforms[0]),
                    transform_cache.lookup(&self.current_transforms[1]),
                ];
                let animated_object2world = AnimatedTransform::new(
                    Arc::clone(&obj2world[0]),
                    Arc::clone(&obj2world[1]),
                    self.render_options.transform_start_time,
                    self.render_options.transform_end_time,
                );
                if prims.len() > 1 {
                    let bvh = BVHAccel::new(&prims, 1, SplitMethod::SAH);
                    prims = vec![Arc::new(bvh)];
                }
                if prims.len() == 1 {
                    let prim = Arc::new(TransformedPrimitive::new(
                        Arc::clone(&prims[0]),
                        animated_object2world,
                    ));
                    prims[0] = prim;
                } else {
                    error!("Error creating TransformedPrimitive in pbrt_shape.");
                }
            }

            // Add `prims` and `area_lights` to scene or current instance.
            if let Some(mut current_instance) = self.render_options.current_instance.clone() {
                if !area_lights.is_empty() {
                    warn!("Area lights not supported with object instancing.");
                    Arc::get_mut(&mut current_instance)
                        .unwrap()
                        .append(&mut prims);
                }
                Arc::get_mut(&mut current_instance)
                    .unwrap()
                    .append(&mut prims);
            } else {
                self.render_options.primitives.append(&mut prims);
                if !area_lights.is_empty() {
                    self.render_options.lights.append(&mut area_lights);
                }
            }
        }
    }

    /// Reverse the orientation of surface normals for shapes that follow this
    /// directive.
    pub fn pbrt_reverse_orientation(&mut self) {
        if self.verify_world("ReverseOrientation") {
            self.graphics_state.reverse_orientation = !self.graphics_state.reverse_orientation;
        }
    }

    /// Begin the definition of a named object instance.
    ///
    /// * `name` - The object instance name.
    pub fn pbrt_object_begin(&mut self, name: String) {
        if self.verify_world("ObjectBegin") {
            self.pbrt_attribute_begin();

            if let Some(_current_instance) = self.render_options.current_instance.clone() {
                error!("ObjectBegin called inside of an instance definition.");
            } else {
                let new_instance: Arc<Vec<ArcPrimitive>> = Arc::new(vec![]);
                self.render_options
                    .instances
                    .insert(name, Arc::clone(&new_instance));
                self.render_options.current_instance = Some(Arc::clone(&new_instance));
            }
        }
    }

    /// End the definition of a named object instance.
    pub fn pbrt_object_end(&mut self) {
        if self.verify_world("ObjectEnd") {
            if let Some(_current_instance) = self.render_options.current_instance.clone() {
                error!("ObjectEnd called outside of instance definition.");
            }
            self.render_options.current_instance = None;

            self.pbrt_attribute_end();
        }
    }

    /// Instantiate a named object.
    ///
    /// * `name` - The object instance name.
    pub fn pbrt_object_instance(&mut self, name: String) {
        if self.verify_world("ObjectInstance") {
            // Perform object instance error checking.
            if let Some(_current_instance) = self.render_options.current_instance.clone() {
                error!("ObjectInstance can't be called inside of instance definition.");
                return;
            }
            if let Some(instance) = self.render_options.instances.get(&name).cloned() {
                let inst = match instance.len() {
                    0 => {
                        return;
                    }
                    1 => Arc::clone(&(&*instance)[0]),
                    _ => {
                        // Create an aggregate for the instance `Primitives`.
                        match GraphicsState::make_accelerator(
                            &self.render_options.accelerator_name,
                            &*instance,
                            &self.render_options.accelerator_params,
                        ) {
                            Ok(acc) => acc.clone(),
                            Err(err) => {
                                error!("{}", err);
                                return;
                            }
                        }
                    }
                };

                // Create `animated_instance_to_world` transform for instance.
                let mut transform_cache = self.transform_cache.lock().unwrap();
                let instance2world = [
                    transform_cache.lookup(&self.current_transforms[0]),
                    transform_cache.lookup(&self.current_transforms[1]),
                ];
                let animated_instance2world = AnimatedTransform::new(
                    Arc::clone(&instance2world[0]),
                    Arc::clone(&instance2world[1]),
                    self.render_options.transform_start_time,
                    self.render_options.transform_end_time,
                );
                let prim = TransformedPrimitive::new(inst, animated_instance2world);
                self.render_options.primitives.push(Arc::new(prim));
            } else {
                error!("Unable to find object instance named '{}'", name);
            }
        }
    }

    /* Helpers */

    /// Returns `true` if the API state is initialized; otherwise it reports
    /// an error and returns `false`.
    ///
    /// * `func` - Function name to report.
    fn verify_initialized(&self, func: &str) -> bool {
        if self.current_api_state == ApiState::Uninitialized {
            error!("pbrt_init() must be before calling '{}()'. Ignoring.", func);
            false
        } else {
            true
        }
    }

    /// Returns `true` if the API state is initialized and inside an Options
    /// block; otherwise it reports an error and returns `false`.
    ///
    /// * `func` - Function name to report.
    fn verify_options(&self, func: &str) -> bool {
        if !self.verify_initialized(func) {
            false
        } else if self.current_api_state == ApiState::WorldBlock {
            error!(
                "Options cannot be set inside world block; 
                '{}()' not allowed. Ignoring.",
                func
            );
            false
        } else {
            true
        }
    }

    /// Returns `true` if the API state is initialized and not inside an Options
    /// block; otherwise it reports an error and returns `false`.
    ///
    /// * `func` - Function name to report.
    fn verify_world(&self, func: &str) -> bool {
        if !self.verify_initialized(func) {
            false
        } else if self.current_api_state == ApiState::OptionsBlock {
            error!(
                "Scene description must beinside world block; 
                '{}()' not allowed. Ignoring.",
                func
            );
            false
        } else {
            true
        }
    }

    /// Emits a warning message if current transformation matrix is for animation.
    ///
    /// * `func` - Function name to report.
    fn warn_if_animated_transform(&self, func: &str) {
        if self.current_transforms.is_animated() {
            warn!(
                "Animated transformations set; ignoring for '{}' 
                and using the start transform only",
                func
            );
        }
    }

    /// Returns a named medium or `None` for the given name.
    ///
    /// * `name` - Medium name.
    /// * `side` - Used to report an error if medium not found.
    fn get_named_medium(&self, name: Option<String>, side: &str) -> Option<ArcMedium> {
        match name {
            Some(n) => {
                if n.is_empty() {
                    error!("Medium name is empty string for side '{}'.", side);
                    None
                } else if let Some(medium) = self.render_options.named_media.get(&n) {
                    Some(medium.clone())
                } else {
                    error!("Named medium '{}' undefined for side '{}'.", n, side);
                    None
                }
            }
            None => None,
        }
    }

    /// Creates a new medium interface.
    fn create_medium_interface(&self) -> MediumInterface {
        let inside =
            self.get_named_medium(self.graphics_state.current_inside_medium.clone(), "inside");
        let outside = self.get_named_medium(
            self.graphics_state.current_outside_medium.clone(),
            "outside",
        );
        MediumInterface::new(inside, outside)
    }
}
