//! Bounding Volume Hierarchy.

use core::geometry::*;
use core::interaction::*;
use core::light::*;
use core::material::*;
use core::paramset::*;
use core::primitive::*;
use core::reflection::*;
use core::stat_inc;

mod common;
mod hlbvh;
mod morton;
mod sah;

pub use common::*;
use shared_arena::{ArenaArc, SharedArena};
use std::sync::{Arc, Mutex};

/// Bounding Volume Hierarchy Accelerator.
#[derive(Clone)]
pub struct BVHAccel {
    /// The primitives in the node.
    pub primitives: Vec<ArcPrimitive>,

    /// Maximum number of primitives in the node. **NOTE**: `u8` limits maximum number to 255.
    pub max_prims_in_node: u8,

    /// Spliting method.
    pub split_method: SplitMethod,

    /// The list of nodes.
    pub nodes: Vec<LinearBVHNode>,
}

impl BVHAccel {
    /// Create a new Bounding Volume Hierarchy Accelerator.
    ///
    /// * `primitives`        - The primitives.
    /// * `max_prims_in_node` - Maximum number of primitives in a node.
    /// * `split_method`      - The splitting method.
    pub fn new(primitives: &[ArcPrimitive], max_prims_in_node: u8, split_method: SplitMethod) -> Self {
        register_stats();

        let n_primitives = primitives.len();
        if n_primitives == 0 {
            Self {
                primitives: primitives.to_vec(),
                max_prims_in_node,
                split_method,
                nodes: vec![],
            }
        } else {
            // Build BVH from primitives.

            // Initializes primitive_info array for primitives.
            let mut primitive_info: Vec<BVHPrimitiveInfo> = Vec::with_capacity(n_primitives);
            for (i, p) in primitives.iter().enumerate() {
                primitive_info.push(BVHPrimitiveInfo::new(i, p.world_bound()));
            }

            // Build BVH tree for primitives using `primitive_info`.
            let arena = SharedArena::<BVHBuildNode>::with_capacity(1024 * 1024);
            let mut total_nodes = 0;
            let ordered_prims = Arc::new(Mutex::new(Vec::<ArcPrimitive>::with_capacity(n_primitives)));

            let root = match split_method {
                SplitMethod::HLBVH => hlbvh::build(
                    &arena,
                    primitives,
                    max_prims_in_node,
                    &mut primitive_info,
                    &mut total_nodes,
                    Arc::clone(&ordered_prims),
                ),
                _ => sah::build(
                    &arena,
                    primitives,
                    split_method,
                    max_prims_in_node,
                    &mut primitive_info,
                    0,
                    n_primitives,
                    &mut total_nodes,
                    Arc::clone(&ordered_prims),
                ),
            };

            let (arena_used, _arean_free) = arena.stats();
            info!(
                "BVH created with {} nodes for {} primitives ({:2} MB), arena allocated {:2} MB",
                total_nodes,
                primitives.len(),
                (total_nodes * std::mem::size_of::<LinearBVHNode>()) as f32 / (1024.0 * 1024.0),
                arena_used as f32 / (1024.0 * 1024.0)
            );

            // Compute representation of depth-first traversal of BVH tree.
            let tree_bytes = total_nodes * std::mem::size_of::<LinearBVHNode>()
                + std::mem::size_of::<Self>()
                + primitives.len() * std::mem::size_of::<ArcPrimitive>();
            stat_inc!(TREE_BYTES, tree_bytes as u64);

            let mut nodes = vec![LinearBVHNode::default(); total_nodes];
            let mut offset = 0_u32;
            Self::flatten_bvh_tree(root, &mut nodes, &mut offset);

            debug_assert!(total_nodes == offset as usize);

            let prims = Arc::clone(&ordered_prims);
            let prims2 = prims.lock().expect("unabled to lock ordered_prims");
            BVHAccel {
                primitives: prims2.to_vec(),
                max_prims_in_node,
                split_method,
                nodes,
            }
        }
    }

    /// Flatten the tree to the linear representation.
    ///
    /// * `node`   - The node.
    /// * `offset` - Tracks current offset into `BVHAccel::nodes`.
    fn flatten_bvh_tree(node: ArenaArc<BVHBuildNode>, nodes: &mut Vec<LinearBVHNode>, offset: &mut u32) -> u32 {
        let my_offset = *offset;
        *offset += 1;

        if node.n_primitives > 0 {
            // This assert!() doesn't make sense. When you have a single primitive in the BVH, children node will be
            // None.
            // assert!(!node.children[0].is_none() && !node.children[1].is_none());
            assert!(node.n_primitives < 65536);

            nodes[my_offset as usize] =
                LinearBVHNode::new_leaf_node(node.bounds, node.first_prim_offset as u32, node.n_primitives as u16);
        } else {
            // Create interior flattened BVH nodes.
            if let Some(child) = node.children[0].clone() {
                // Ignore first child offset for interior node.
                Self::flatten_bvh_tree(child, nodes, offset);
            }

            if let Some(child) = node.children[1].clone() {
                let second_child_offset = Self::flatten_bvh_tree(child, nodes, offset);
                nodes[my_offset as usize] =
                    LinearBVHNode::new_interior_node(node.bounds, second_child_offset as u32, node.split_axis.into());
            }
        }

        my_offset
    }
}

/// Tag `BVHAccel` as an `Aggregate`.
impl Aggregate for BVHAccel {}

impl Primitive for BVHAccel {
    /// Returns a bounding box in the world space.
    fn world_bound(&self) -> Bounds3f {
        if !self.nodes.is_empty() {
            self.nodes[0].bounds
        } else {
            Bounds3f::EMPTY
        }
    }

    /// Returns geometric details if a ray intersects the primitive and updates the t_max parameter of the ray. If there
    /// is no intersection, `None` is returned.
    ///
    /// * `r` - The ray.
    fn intersect(&self, r: &mut Ray) -> Option<SurfaceInteraction> {
        let mut si: Option<SurfaceInteraction> = None;
        if !self.nodes.is_empty() {
            let inv_dir = Vector3f::new(1.0 / r.d.x, 1.0 / r.d.y, 1.0 / r.d.z);
            let dir_is_neg = [
                if inv_dir.x < 0.0 { 1_u8 } else { 0_u8 },
                if inv_dir.y < 0.0 { 1_u8 } else { 0_u8 },
                if inv_dir.z < 0.0 { 1_u8 } else { 0_u8 },
            ];

            // Follow ray through BVH nodes to find primitive intersections.
            let (mut to_visit_offset, mut current_node_index) = (0, 0);
            let mut nodes_to_visit = [0_usize; 64];

            loop {
                // Check ray against BVH node.
                let node = &self.nodes[current_node_index];
                if node.bounds.intersect_p_inv(r, &inv_dir, dir_is_neg) {
                    if node.n_primitives > 0 {
                        // Intersect ray with primitives in leaf BVH node.
                        for i in 0..node.n_primitives {
                            let idx = node.offset as usize + i as usize;
                            if let Some(hit) = self.primitives[idx].intersect(r) {
                                si = Some(hit);
                            }
                        }
                        if to_visit_offset == 0 {
                            break;
                        }
                        to_visit_offset -= 1;
                        current_node_index = nodes_to_visit[to_visit_offset];
                    } else {
                        // Put far BVH node on nodes_to_visit stack, advance to near node.
                        if dir_is_neg[node.axis as usize] == 1 {
                            nodes_to_visit[to_visit_offset] = current_node_index + 1;
                            to_visit_offset += 1;
                            current_node_index = node.offset as usize;
                        } else {
                            nodes_to_visit[to_visit_offset] = node.offset as usize;
                            to_visit_offset += 1;
                            current_node_index += 1;
                        }
                    }
                } else {
                    if to_visit_offset == 0 {
                        break;
                    }
                    to_visit_offset -= 1;
                    current_node_index = nodes_to_visit[to_visit_offset];
                }
            }
        }
        si
    }

    /// Returns `true` if a ray-primitive intersection succeeds; otherwise `false`.
    ///
    /// * `r` - The ray.
    fn intersect_p(&self, r: &Ray) -> bool {
        if !self.nodes.is_empty() {
            let inv_dir = Vector3f::new(1.0 / r.d.x, 1.0 / r.d.y, 1.0 / r.d.z);
            let dir_is_neg = [
                if inv_dir.x < 0.0 { 1_u8 } else { 0_u8 },
                if inv_dir.y < 0.0 { 1_u8 } else { 0_u8 },
                if inv_dir.z < 0.0 { 1_u8 } else { 0_u8 },
            ];

            // Follow ray through BVH nodes to find primitive intersections.
            let (mut to_visit_offset, mut current_node_index) = (0, 0);
            let mut nodes_to_visit = [0_usize; 64];

            loop {
                // Check ray against BVH node
                let node = &self.nodes[current_node_index];
                if node.bounds.intersect_p_inv(r, &inv_dir, dir_is_neg) {
                    if node.n_primitives > 0 {
                        // Intersect ray with primitives in leaf BVH node.
                        for i in 0..node.n_primitives {
                            let idx = node.offset as usize + i as usize;
                            if self.primitives[idx].intersect_p(r) {
                                return true;
                            }
                        }
                        if to_visit_offset == 0 {
                            break;
                        }
                        to_visit_offset -= 1;
                        current_node_index = nodes_to_visit[to_visit_offset];
                    } else {
                        // Put far BVH node on nodes_to_visit stack, advance to near node.
                        if dir_is_neg[node.axis as usize] == 1 {
                            nodes_to_visit[to_visit_offset] = current_node_index + 1;
                            to_visit_offset += 1;
                            current_node_index = node.offset as usize;
                        } else {
                            nodes_to_visit[to_visit_offset] = node.offset as usize;
                            to_visit_offset += 1;
                            current_node_index += 1;
                        }
                    }
                } else {
                    if to_visit_offset == 0 {
                        break;
                    }
                    to_visit_offset -= 1;
                    current_node_index = nodes_to_visit[to_visit_offset];
                }
            }
        }
        false
    }

    /// Returns a reference to the AreaLight that describes the primitive’s emission distribution, if the primitive is
    /// itself a light source. If the primitive is not emissive, this method should return `None`.
    ///
    /// *NOTE*: This should never be called. Calling code should directly call get_area_light() on the primitive from
    /// the ray-primitive intersection.
    fn get_area_light(&self) -> Option<ArcLight> {
        error!(
            "TransformedPrimitive::get_area_light() shouldn't be called; \
            should've gone to GeometricPrimitive."
        );
        None
    }

    /// Returns a reference to the material instance assigned to the primitive. If `None` is returned, ray intersections
    /// with the primitive should be ignored; the primitive only serves to delineate a volume of space for participating
    /// media. This method is also used to check if two rays have intersected the same object by comparing their
    /// Material pointers.
    ///
    /// *NOTE*: This should never be called. Calling code should directly call get_material() on the primitive from
    /// the ray-primitive intersection.
    fn get_material(&self) -> Option<ArcMaterial> {
        error!(
            "TransformedPrimitive::get_material() shouldn't be called; \
            should've gone to GeometricPrimitive."
        );
        None
    }

    /// Initializes representations of the light-scattering properties of the material at the intersection point on the
    /// surface.
    ///
    /// *NOTE*: This should never be called. Calling code should directly call compute_scattering_functions() on the
    /// primitive from the ray-primitive intersection.
    ///
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode.
    /// * `allow_multiple_lobes` - Allow multiple lobes.
    /// * `bsdf`                 - The computed BSDF.
    /// * `bssrdf`               - The computed BSSSRDF.
    fn compute_scattering_functions<'scene>(
        &self,
        _si: &mut SurfaceInteraction<'scene>,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
        _bsdf: &mut Option<BSDF>,
        _bssrdf: &mut Option<BSDF>,
    ) {
        error!(
            "TransformedPrimitive::compute_scattering_functions() shouldn't be \
            called; should've gone to GeometricPrimitive."
        );
    }
}

impl From<(&ParamSet, &[ArcPrimitive])> for BVHAccel {
    /// Create a `BVHAccel` from given parameter set and primitives.
    ///
    /// * `p` - Tuple containing the parameter set and primitives.
    fn from(p: (&ParamSet, &[ArcPrimitive])) -> Self {
        let (params, prims) = p;
        let split_method_name = params.find_one_string("splitmethod", String::from("sah"));
        let split_method = match &split_method_name[..] {
            "sah" => SplitMethod::SAH,
            "hlbvh" => SplitMethod::HLBVH,
            "middle" => SplitMethod::Middle,
            "equal" => SplitMethod::EqualCounts,
            sm => {
                warn!("BVH split method '{}' unknown.  Using 'sah'.", sm);
                SplitMethod::SAH
            }
        };

        let max_prims_in_node = params.find_one_int("maxnodeprims", 4) as u8;
        Self::new(prims, max_prims_in_node, split_method)
    }
}
