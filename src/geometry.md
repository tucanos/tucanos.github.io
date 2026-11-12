# Geometry

Relying on the mesh boundary alone to describe the geometry is often not sufficient to ensure consistency across mesh adaptation loops. The `Geometry` trait is used to:

- validate topology compatibility (`check(topo)`),
- project points onto geometric support (`project(pt, tag)`),
- compare mesh normals to geometry normals (`angle(pt, n, tag)`).

Projection and normal angle calculations are performed for a given tag, not just for nearest surface globally. This is how `tucanos` preserves sharp features and interfaces during adaptation, and this is why the mesh and the geometry must have consistent topology hierarchies.

## Meshed geometries

In `tucanos`, the main representation of geometries is through tesselations as it enables fast projections and can be used together with mesh deformation workflows. `MeshedGeometry` stores one spatial index per tagged patch in order to enable tag-routed projection and normal queries (using the [Element/object index](spatial_indices.md#elementobject-index)).

> [!NOTE]
> In practice, many inputs are piecewise-linear surfaces where triangles can be extremely anisotropic (very flat/slender). The spatial index provided by parry handles this robustly.

Face tags should be consistent between the geometry and the mesh. However edge (in 3D) and node tags are built automatically and it is important to ensure consistency for edge projection. This is why
- edge entities are identified by their parent face-tag set,
- `set_topo_map(...)` maps mesh edge tags to geometry edge tags from those parents,
- projection on edge tags uses this map before querying edge geometry.

This ensures that an edge/node in the adapted mesh is projected onto the correct geometric feature, even if mesh-local tag ids differ from geometry-local ids.

> [!TIP]
> **Geometry tag lookup errors**
>
> `MeshedGeometry::new(...)`, `set_topo_map(...)`, `project(...)`, and `angle(...)` panic when a tag cannot be resolved:
>
> | Message | Meaning / action |
> | --- | --- |
> | `Invalid face tag {tag:?}. Available tags are {:?}` | `project`/`angle` was called with a face tag that is not present in the `MeshedGeometry`. Build the geometry from the same tagged boundary and keep mesh/geometry face tags consistent. |
> | `Invalid edge tag {tag:?}. Available tags are {:?}` | An edge tag could not be mapped to a geometry edge. Check that edge tags are inferred consistently and that `set_topo_map(...)` was called. |
> | `Multiple edge tags with parents {parents:?}` | `set_topo_map` found more than one mesh edge tag for the same parent face-tag set. The mesh and geometry topologies are inconsistent. |
> | `Geometry mesh empty for tag {tag}` | A face tag was found on the boundary mesh but selected no face when building its patch. |
>
> `check(topo)` does not re-validate these mappings (`MeshedGeometry` checks them when it is built); a panic therefore indicates inconsistent input tags or an incomplete topology map.

> [!TIP]
> `MeshedGeometry` can be used with linear or quadratic meshes.

### Curvature handling in `MeshedGeometry`

Curvature is precomputed patch-by-patch (for faces) when `MeshedGeometry` is built. The curvature at a point `pt` for tag `tag` can then be obtained with `geom.curvature(pt, tag)` that performs the following operations:
1. select the patch from `tag`,
2. find the nearest patch element to `pt`,
3. return that element's stored principal curvature direction(s).

Returned values follow the same convention as the curvature module:
- in 2D curves: one vector `u`,
- in 3D surfaces: two vectors `(u, v)`.

Each vector direction gives a principal direction and its norm gives the associated principal curvature magnitude (see [here](curvature.md))

The same per-patch data can be exported with `write_curvature(...)` for inspection/debugging.


### Rust example: 

```rust
# extern crate tmesh;
# extern crate tucanos;
use tmesh::{Tag, Vert3d};
use tmesh::mesh::{BoundaryMesh3d, Mesh, Mesh3d, box_mesh};
use tucanos::geometry::{Geometry, MeshedGeometry};
use tucanos::mesh::MeshTopology;

// Build a tetrahedral box mesh
let msh = box_mesh::<Mesh3d>(1.0, 6, 1.0, 6, 1.0, 6);

// Extract its boundary surface mesh and split it (in practice, one would use a 
// better representation of the geometry that the mesh boundary)
let (bdy, _) = msh.boundary::<BoundaryMesh3d>();
let mut bdy = bdy.split();
bdy.fix().unwrap();

// Build mesh topology and geometry, then connect edge-tag spaces
let topo = MeshTopology::new(&msh);
let mut geom = MeshedGeometry::new(&bdy).unwrap();
geom.set_topo_map(topo.topo());
```

## Analytical geometries

Analytical geometries (see `CylinderGeometry` in `tucanos/src/mesh/test_meshes.rs`) are mainly used for testing and controlled validation cases. The idea is to provide closed-form projection and normal-angle evaluation without relying on a meshed surface index. This is especially useful to:

- test projection and tagging logic in isolation,
- validate remeshing quality against known exact geometry,
- build reproducible benchmarks.

## Outlook: CAD geometry with [truck](https://github.com/ricosjp/truck)

`truck` could be used as an alternative backend for `Geometry` when the input geometry is CAD-based (for example STEP files) rather than a pre-tessellated boundary mesh.

Potential workflow:

1. Import a STEP model and extract B-Rep entities (faces/edges/vertices).
2. Build a CAD-geometry map from topological tags to one or more CAD patches.
3. Implement `project(pt, tag)` by running nearest-point projection on the mapped CAD patch set.
4. Implement `angle(pt, n, tag)` by evaluating CAD surface normals at projected points.

This would make it possible to project directly on exact CAD surfaces instead of a triangle approximation.

Main practical difficulty:

- matching mesh face tags to CAD patches is not always one-to-one;
- in real CAD models, one mesh tag can correspond to several CAD patches (trimmed surfaces, split faces, imported topology differences), so projection must search across a patch set per tag.

Performance caveat:

- CAD projection is computationally more expensive than projection on a triangle mesh with a prebuilt object index;
- in adaptation loops with many vertex projections, this cost can become a dominant runtime term unless caching, locality-based initialization, or hybrid strategies are used.


Illustrative `truck` snippets (marked `ignore`) are shown below.

### STEP import with `truck`

```rust,ignore
use std::fs;

// Example only: exact modules/types may differ between truck versions.
// Goal: read STEP text, parse it, and obtain a topological shell/solid.
let step_text = fs::read_to_string("model.step")?;
let step_model = truck_stepio::r#in::step_to_model(&step_text)?;

// Extract CAD topology/geometry container used by the downstream code.
let cad_shape = step_model.to_compressed_shell()?;
```

### Point projection on CAD patches

```rust,ignore
use truck_base::cgmath::Point3;

// Example only: tag -> candidate CAD patches map is precomputed.
let p = Point3::new(0.2, 0.3, 0.7);
let candidate_patches = tag_to_patches.get(&face_tag).unwrap();

// Project on all candidate patches and keep the nearest result.
let mut best = None;
for patch in candidate_patches {
	// Returns (u, v, projected_point, distance)
	let hit = patch.search_nearest_parameter(p)?;
	best = match best {
		None => Some(hit),
		Some(prev) if hit.distance < prev.distance => Some(hit),
		_ => best,
	};
}

let projected = best.unwrap().projected_point;
```

### Normal computation at projected point

```rust,ignore
use truck_base::cgmath::InnerSpace;

// Suppose projection returned (patch_id, u, v, projected_point).
let patch = &cad_patches[patch_id];

// Evaluate CAD surface derivatives and compute a unit normal.
let du = patch.uder(u, v);
let dv = patch.vder(u, v);
let n_cad = du.cross(dv).normalize();

// Compare with mesh normal for angle-based quality control.
let cos_a = mesh_normal.dot(n_cad).clamp(-1.0, 1.0);
let angle_deg = cos_a.acos().to_degrees();
```