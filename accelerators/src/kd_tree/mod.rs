//! KD Tree.

use core::geometry::*;
use core::light::*;
use core::material::*;
use core::paramset::*;
use core::pbrt::*;
use core::primitive::*;
use std::sync::Arc;

mod common;
use common::*;

/// KD Tree Accelerator.
#[derive(Clone)]
pub struct KDTreeAccel {
    /// Intersection cost.
    pub isect_cost: i32,

    /// Traversal cost.
    pub traversal_cost: i32,

    /// Maximum number of primitives in leaf node.
    pub max_prims: u32,

    /// Bonus value used when one of the 2 regions along split are empty.
    pub empty_bonus: Float,

    /// The primitives.
    pub primitives: Vec<ArcPrimitive>,

    /// Indexes into `primitives`.
    pub primitive_indices: Vec<u32>,

    /// The tree root.
    pub nodes: Vec<KdAccelNode>,

    /// Number of allocated nodes.
    pub n_alloced_nodes: u32,

    /// Next free node.
    pub next_free_node: u32,

    /// Bounding box.
    pub bounds: Bounds3f,
}

impl KDTreeAccel {
    /// Create a new KD Tree Accelerator.
    ///
    /// * `primitives`     - Pimitives.
    /// * `isect_cost`     - Intersection cost.
    /// * `traversal_cost` - Traversal cost.
    /// * `empty_bonus`    - Bonus value used when one of the 2 regions along
    ///                      split are empty.
    /// * `max_prims`      - Maximum number of primitives in leaf.
    /// * `max_depth`      - Maximum depth of tree.
    pub fn new(
        primitives: &[ArcPrimitive],
        isect_cost: i32,
        traversal_cost: i32,
        empty_bonus: Float,
        max_prims: u32,
        max_depth: i32,
    ) -> Self {
        // Build kd-tree for accelerator.
        let count = primitives.len();
        let next_free_node = 0;
        let n_alloced_nodes = 0;

        let max_depth = if max_depth < 0 {
            (8.0 + 1.3 * Log2::log2(count as i64) as Float).round() as i32
        } else {
            max_depth
        };

        // Compute bounds for kd-tree construction.
        let mut bounds = Bounds3f::empty();
        let mut prim_bounds = Vec::<Bounds3f>::with_capacity(count);
        for prim in primitives.iter() {
            let b = prim.world_bound();
            bounds = bounds.union(&b);
            prim_bounds.push(b);
        }

        // Allocate working memory for kd-tree construction.
        let mut edges: [Vec<BoundEdge>; 3] = [
            Vec::<BoundEdge>::with_capacity(2 * count),
            Vec::<BoundEdge>::with_capacity(2 * count),
            Vec::<BoundEdge>::with_capacity(2 * count),
        ];

        let mut prims0 = Vec::<u32>::with_capacity(count);
        let mut prims1: Vec<u32> = vec![0_u32; (max_depth + 1) as usize * count];

        // Initialize prim_nums for kd-tree construction.
        let prim_nums: Vec<u32> = (0..count as u32).collect();

        // Start recursive construction of kd-tree
        let mut kd_tree = Self {
            isect_cost,
            traversal_cost,
            max_prims,
            empty_bonus,
            primitives: primitives.to_vec(),
            primitive_indices: vec![],
            nodes: vec![],
            n_alloced_nodes,
            next_free_node,
            bounds,
        };

        kd_tree.build_tree(
            0,
            &bounds,
            &prim_bounds,
            &prim_nums[..],
            count as u32,
            max_depth,
            &mut edges,
            &mut prims0[..],
            &mut prims1[..],
            0,
        );

        kd_tree
    }

    fn build_tree(
        &mut self,
        node_num: u32,
        node_bounds: &Bounds3f,
        all_prim_bounds: &[Bounds3f],
        prim_nums: &[u32],
        n_primitives: u32,
        depth: i32,
        edges: &mut [Vec<BoundEdge>; 3],
        prims0: &mut [u32],
        prims1: &mut [u32],
        bad_refines: i32,
    ) {
        debug_assert!(node_num == self.next_free_node);

        // Get next free node from nodes.
        if self.next_free_node == self.n_alloced_nodes {
            let n_new_alloc_nodes = max(2 * self.n_alloced_nodes, 512);
            self.nodes
                .resize_with(n_new_alloc_nodes as usize, KdAccelNode::default);
            self.n_alloced_nodes = n_new_alloc_nodes;
        }
        self.next_free_node += 1;

        // Initialize leaf node if termination criteria met.
        if n_primitives <= self.max_prims || depth == 0 {
            self.nodes[node_num as usize].init_leaf(
                prim_nums,
                n_primitives,
                &mut self.primitive_indices,
            );
        }

        // Initialize interior node and continue recursion.

        // Choose split axis position for interior node.
        let mut best_axis: Option<usize> = None;
        let mut best_offset: Option<usize> = None;
        let mut best_cost = INFINITY;
        let old_cost = self.isect_cost as Float * n_primitives as Float;
        let total_sa = node_bounds.surface_area();
        let inv_total_sa = 1.0 / total_sa;
        let d = node_bounds.diagonal();

        // Choose which axis to split.
        let mut axis: usize = node_bounds.maximum_extent().into();
        let mut retries = 0;

        loop {
            // Initialize edges for `axis`.
            for (i, pn) in prim_nums.iter().enumerate().take(n_primitives as usize) {
                let bounds = all_prim_bounds[*pn as usize];
                edges[axis][2 * i] = BoundEdge::new(bounds.p_min[axis], *pn, true);
                edges[axis][2 * i + 1] = BoundEdge::new(bounds.p_max[axis], *pn, false);
            }

            // Sort `edges` for `axis`.
            edges[axis].sort_by(|e0, e1| {
                if e0.t == e1.t {
                    e0.edge_type.partial_cmp(&e1.edge_type).unwrap()
                } else {
                    e0.t.partial_cmp(&e1.t).unwrap()
                }
            });

            // Compute cost of all splits for `axis` to find best.
            let (mut n_below, mut n_above) = (0, n_primitives);
            for i in 0..2 * n_primitives as usize {
                if edges[axis][i].edge_type == EdgeType::End {
                    n_above -= 1;
                }

                let edge_t = edges[axis][i].t;
                if edge_t > node_bounds.p_min[axis] && edge_t < node_bounds.p_max[axis] {
                    // Compute cost for split at i^th edge.

                    // Compute child surface areas for split at `edge_t`.
                    let other_axis_0 = (axis + 1) % 3;
                    let other_axis_1 = (axis + 2) % 3;
                    let below_sa = 2.0
                        * (d[other_axis_0] * d[other_axis_1]
                            + (edge_t - node_bounds.p_min[axis])
                                * (d[other_axis_0] + d[other_axis_1]));
                    let above_sa = 2.0
                        * (d[other_axis_0] * d[other_axis_1]
                            + (node_bounds.p_max[axis] - edge_t)
                                * (d[other_axis_0] + d[other_axis_1]));
                    let p_below = below_sa * inv_total_sa;
                    let p_above = above_sa * inv_total_sa;
                    let eb = if n_above == 0 || n_below == 0 {
                        self.empty_bonus
                    } else {
                        0.0
                    };
                    let cost = self.traversal_cost as Float
                        + (self.isect_cost as Float)
                            * (1.0 - eb)
                            * (p_below * n_below as Float + p_above * n_above as Float);

                    // Update best split if this is lowest cost so far
                    if cost < best_cost {
                        best_cost = cost;
                        best_axis = Some(axis);
                        best_offset = Some(i);
                    }
                }
                if edges[axis][i].edge_type == EdgeType::Start {
                    n_below += 1;
                }
            }
            debug_assert!(n_below == n_primitives && n_above == 0);

            // Create leaf node if no good splits were found.
            if best_axis.is_none() && retries < 2 {
                retries += 1;
                axis = (axis + 1) % 3;
                continue;
            }

            let new_bad_refines = if best_cost > old_cost {
                bad_refines + 1
            } else {
                bad_refines
            };

            if (best_cost > 4.0 * old_cost && n_primitives < 16)
                || best_axis.is_none()
                || new_bad_refines == 3
            {
                self.nodes[node_num as usize].init_leaf(
                    prim_nums,
                    n_primitives,
                    &mut self.primitive_indices,
                );
                return;
            }

            // Classify primitives with respect to split
            let (mut n0, mut n1) = (0_u32, 0_u32);
            let best_axis = best_axis.unwrap();
            let best_offset = best_offset.unwrap();

            for i in 0..best_offset {
                if edges[best_axis][i].edge_type == EdgeType::Start {
                    prims0[n0 as usize] = edges[best_axis][i].prim_num;
                    n0 += 1;
                }
            }

            for i in best_offset + 1..2 * n_primitives as usize {
                if edges[best_axis][i].edge_type == EdgeType::End {
                    prims1[n1 as usize] = edges[best_axis][i].prim_num;
                    n1 += 1;
                }
            }

            // Recursively initialize children nodes
            let t_split = edges[best_axis][best_offset].t;

            let mut bounds0 = *node_bounds;
            let mut bounds1 = *node_bounds;
            bounds0.p_max[best_axis] = t_split;
            bounds1.p_min[best_axis] = t_split;

            let n = n_primitives as usize;
            self.build_tree(
                node_num + 1,
                &bounds0,
                all_prim_bounds,
                prim_nums,
                n0,
                depth - 1,
                edges,
                prims0,
                &mut prims1[n..],
                new_bad_refines,
            );

            let above_child = self.next_free_node;
            self.nodes[node_num as usize].init_interior(best_axis as u32, above_child, t_split);

            self.build_tree(
                above_child,
                &bounds1,
                all_prim_bounds,
                prim_nums,
                n1,
                depth - 1,
                edges,
                prims0,
                &mut prims1[n..],
                new_bad_refines,
            );
        }
    }
}

/// Tag `KDTreeAccel` as an `Aggregate`.
impl Aggregate for KDTreeAccel {}

impl Primitive for KDTreeAccel {
    /// Returns a bounding box in the world space.
    fn world_bound(&self) -> Bounds3f {
        self.bounds
    }

    /// Returns geometric details if a ray intersects the primitive and updates
    /// the t_max parameter of the ray. If there is no intersection, `None` is
    /// returned.
    ///
    /// * `r`                  - The ray.
    fn intersect(&self, r: &mut Ray) -> Option<SurfaceInteraction> {
        let mut si: Option<SurfaceInteraction> = None;

        // Compute initial parametric range of ray inside kd-tree extent.
        if let Some((mut t_min, mut t_max)) = self.bounds.intersect_p(r) {
            // Prepaer to traverse kd-tree for ray.
            let inv_dir = Vector3f::new(1.0 / r.d.x, 1.0 / r.d.y, 1.0 / r.d.z);
            let mut todo = [KdToDo::default(); MAX_TO_DO];
            let mut todo_pos = 0;
            let mut current_node_idx = 0;
            let mut current_node = self.nodes.get(current_node_idx);

            // Traverse kd-tree nodes in order for ray.
            while let Some(node) = current_node {
                // Bail out if we found a hit closer than the current node.
                if r.t_max < t_min {
                    break;
                }

                if !node.is_leaf() {
                    // Process kd-tree interior node.

                    // Compute parametric distance along ray to split plane.
                    let axis = node.split_axis() as usize;
                    let t_plane = (node.split_pos() - r.o[axis]) * inv_dir[axis];

                    // Get node children pointers for ray.
                    let below_first = (r.o[axis] < node.split_pos())
                        || (r.o[axis] == node.split_pos() && r.d[axis] <= 0.0);

                    let (first_child_idx, second_child_idx) = if below_first {
                        (current_node_idx + 1, node.above_child() as usize)
                    } else {
                        (node.above_child() as usize, current_node_idx + 1)
                    };
                    let first_child = self.nodes.get(first_child_idx);
                    let second_child = self.nodes.get(second_child_idx);

                    // Advance to next child node, possibly enqueue other child.
                    if t_plane > t_max || t_plane <= 0.0 {
                        current_node = first_child;
                        current_node_idx = first_child_idx;
                    } else if t_plane < t_min {
                        current_node = second_child;
                        current_node_idx = second_child_idx;
                    } else {
                        // Enqueue _second_child_ in todo list.
                        todo[todo_pos].node = second_child;
                        todo[todo_pos].node_idx = second_child_idx;
                        todo[todo_pos].t_min = t_plane;
                        todo[todo_pos].t_max = t_max;
                        todo_pos += 1;

                        current_node = first_child;
                        current_node_idx = first_child_idx;
                        t_max = t_plane;
                    }
                } else {
                    // Check for intersections inside leaf node.
                    let n_primitives = node.n_primitives();
                    if n_primitives == 1 {
                        // Check one primitive inside leaf node.
                        let one_primitive = node.one_primitive() as usize;
                        si = self.primitives[one_primitive].intersect(r);
                    } else {
                        for i in 0..n_primitives as usize {
                            // Check one primitive inside leaf node.
                            let offset = node.primitive_indices_offset() as usize;
                            let index = self.primitive_indices[offset + i] as usize;
                            si = self.primitives[index].intersect(r);
                        }
                    }

                    // Grab next node to process from todo list
                    if todo_pos > 0 {
                        todo_pos -= 1;
                        current_node = todo[todo_pos].node;
                        current_node_idx = todo[todo_pos].node_idx;
                        t_min = todo[todo_pos].t_min;
                        t_max = todo[todo_pos].t_max;
                    } else {
                        break;
                    }
                }
            }
        }

        si
    }

    /// Returns `true` if a ray-primitive intersection succeeds; otherwise `false`.
    ///
    /// * `r`                  - The ray.
    fn intersect_p(&self, r: &Ray) -> bool {
        // Compute initial parametric range of ray inside kd-tree extent.
        if let Some((mut t_min, mut t_max)) = self.bounds.intersect_p(r) {
            // Prepaer to traverse kd-tree for ray.
            let inv_dir = Vector3f::new(1.0 / r.d.x, 1.0 / r.d.y, 1.0 / r.d.z);
            let mut todo = [KdToDo::default(); MAX_TO_DO];
            let mut todo_pos = 0;
            let mut current_node_idx = 0;
            let mut current_node = self.nodes.get(current_node_idx);

            // Traverse kd-tree nodes in order for ray.
            while let Some(node) = current_node {
                // Bail out if we found a hit closer than the current node.
                if r.t_max < t_min {
                    break;
                }

                if !node.is_leaf() {
                    // Process kd-tree interior node.
                    // Compute parametric distance along ray to split plane.
                    let axis = node.split_axis() as usize;
                    let t_plane = (node.split_pos() - r.o[axis]) * inv_dir[axis];

                    // Get node children pointers for ray.
                    let below_first = (r.o[axis] < node.split_pos())
                        || (r.o[axis] == node.split_pos() && r.d[axis] <= 0.0);

                    let (first_child_idx, second_child_idx) = if below_first {
                        (current_node_idx + 1, node.above_child() as usize)
                    } else {
                        (node.above_child() as usize, current_node_idx + 1)
                    };
                    let first_child = self.nodes.get(first_child_idx);
                    let second_child = self.nodes.get(second_child_idx);

                    // Advance to next child node, possibly enqueue other child.
                    if t_plane > t_max || t_plane <= 0.0 {
                        current_node = first_child;
                        current_node_idx = first_child_idx;
                    } else if t_plane < t_min {
                        current_node = second_child;
                        current_node_idx = second_child_idx;
                    } else {
                        // Enqueue `second_child` in todo list.
                        todo[todo_pos].node = second_child;
                        todo[todo_pos].node_idx = second_child_idx;
                        todo[todo_pos].t_min = t_plane;
                        todo[todo_pos].t_max = t_max;
                        todo_pos += 1;

                        current_node = first_child;
                        current_node_idx = first_child_idx;
                        t_max = t_plane;
                    }
                } else {
                    // Check for intersections inside leaf node.
                    let n_primitives = node.n_primitives();
                    if n_primitives == 1 {
                        let one_primitive = node.one_primitive() as usize;
                        let p = Arc::clone(&self.primitives[one_primitive]);

                        // Check one primitive inside leaf node.
                        if p.intersect_p(r) {
                            return true;
                        }
                    } else {
                        for i in 0..n_primitives as usize {
                            let offset = node.primitive_indices_offset() as usize;
                            let index = self.primitive_indices[offset + i] as usize;
                            let p = Arc::clone(&self.primitives[index]);

                            // Check one primitive inside leaf node.
                            if p.intersect_p(r) {
                                return true;
                            }
                        }
                    }

                    // Grab next node to process from todo list
                    if todo_pos > 0 {
                        todo_pos -= 1;
                        current_node = todo[todo_pos].node;
                        current_node_idx = todo[todo_pos].node_idx;
                        t_min = todo[todo_pos].t_min;
                        t_max = todo[todo_pos].t_max;
                    } else {
                        break;
                    }
                }
            }
        }

        false
    }

    /// Returns a reference to the AreaLight that describes the primitive’s
    /// emission distribution, if the primitive is itself a light source.
    /// If the primitive is not emissive, this method should return `None`.
    ///
    /// *NOTE*: This should never be called. Calling code should directly call
    /// get_area_light() on the primitive from the ray-primitive intersection.
    fn get_area_light(&self) -> Option<ArcAreaLight> {
        error!(
            "TransformedPrimitive::get_area_light() shouldn't be called; \
            should've gone to GeometricPrimitive."
        );
        None
    }

    /// Returns a reference to the material instance assigned to the primitive.
    /// If `None` is returned, ray intersections with the primitive should be
    /// ignored; the primitive only serves to delineate a volume of space for
    /// participating media. This method is also used to check if two rays have
    /// intersected the same object by comparing their Material pointers.
    ///
    /// *NOTE*: This should never be called. Calling code should directly call
    /// get_material() on the primitive from the ray-primitive intersection.
    fn get_material(&self) -> Option<ArcMaterial> {
        error!(
            "TransformedPrimitive::get_material() shouldn't be called; \
            should've gone to GeometricPrimitive."
        );
        None
    }

    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// *NOTE*: This should never be called. Calling code should directly call
    /// compute_scattering_functions() on the primitive from the ray-primitive
    /// intersection.
    ///
    /// * `_si`                   - The surface interaction at the intersection.
    /// * `_mode`                 - Transport mode.
    /// * `_allow_multiple_lobes` - Allow multiple lobes.
    fn compute_scattering_functions(
        &self,
        _si: &mut SurfaceInteraction,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        error!(
            "TransformedPrimitive::compute_scattering_functions() shouldn't be \
            called; should've gone to GeometricPrimitive."
        );
    }
}

impl From<(&ParamSet, &[ArcPrimitive])> for KDTreeAccel {
    /// Create a `KDTreeAccel ` from given parameter set and primitives.
    ///
    /// * `p` - Tuple containing the parameter set and primitives.
    fn from(p: (&ParamSet, &[ArcPrimitive])) -> Self {
        let (params, prims) = p;
        let isect_cost = params.find_one_int("intersectcost", 80);
        let trav_cost = params.find_one_int("traversalcost", 1);
        let empty_bonus = params.find_one_float("emptybonus", 0.5);
        let max_prims = params.find_one_int("maxprims", 1) as u32;
        let max_depth = params.find_one_int("maxdepth", -1);

        Self::new(
            prims,
            isect_cost,
            trav_cost,
            empty_bonus,
            max_prims,
            max_depth,
        )
    }
}
