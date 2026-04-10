# Geometry

This module provides the structural foundation for associating computational meshes with their underlying geometric domains. It facilitates projecting mesh vertices onto true geometric surfaces (e.g., during mesh adaptation or order elevation), computing geometric deviations (like distance and normal angles), and handling surface curvature.

---

## Core Trait: `Geometry<const D: usize>`

The `Geometry` trait defines the essential interface for any geometric representation in a `D`-dimensional space. It is designed to be thread-safe (`Send + Sync`).

### Key Methods

* **`check(&self, topo: &Topology) -> Result<()>`**
  Validates that the geometric model is consistent with the provided topological hierarchy.
* **`project(&self, pt: &mut Vertex<D>, tag: &TopoTag) -> f64`**
  Projects a given vertex `pt` onto the geometric entity identified by `tag`. It modifies the vertex coordinates in-place to the projected coordinates and returns the projection distance.
* **`angle(&self, pt: &Vertex<D>, n: &Vertex<D>, tag: &TopoTag) -> f64`**
  Computes the angle (in degrees) between a given normal vector `n` and the true geometric normal at the projection of `pt` on the specified geometric surface.
* **`project_vertices<M>(&self, mesh: &mut M, topo: &MeshTopology) -> f64`**
  Iterates through all vertices in a mesh and projects boundary vertices (entities with a dimension `< D`) onto the geometry, returning the maximum projection distance.
* **`max_distance<M>(&self, mesh: &M) -> f64`**
  Computes the maximum geometric deviation (distance) between the mesh's face/element centers and the true geometry.
* **`max_normal_angle<M>(&self, mesh: &M) -> f64`**
  Computes the maximum angle deviation between the mesh's discrete normals and the true continuous geometry normals.
* **`to_quadratic_triangle_mesh` & `to_quadratic_edge_mesh`**
  Elevates a linear mesh (straight edges/flat faces) to a quadratic mesh. It achieves this by inserting mid-edge nodes and projecting those new nodes onto the true geometry using `project_vertices`, thereby curving the mesh to match the physical boundaries.

---

## Implementations

### 1. `NoGeometry`
A dummy implementation representing the absence of an underlying geometric model. 
* `project()` leaves the vertex unchanged and returns `0.0`.
* `angle()` always returns `0.0`.

### 2. `MeshedGeometry<const D: usize, M: Mesh<D>>`
A discrete, STL-like representation of a geometry constructed from a high-resolution boundary mesh. It relies on spatial indices (`ObjectIndex`, a bounding volume hierarchy) to efficiently find the closest projection points on the discrete surface.

#### Internal Patch Structures
To manage different topological dimensions, the geometry is split into discrete patches:
* **`MeshedPatchGeometry`**: Represents a specific topological edge/curve (dimension `D-1`). It holds the elements belonging to a specific topological tag and uses a spatial index to compute projections.
* **`MeshedPatchGeometryWithCurvature`**: Represents a specific topological surface patch (dimension `D`). In addition to projecting points, it pre-computes and stores principal curvature directions (`u` and optionally `v` for 3D flows) for anisotropic refinement.

#### Attributes of `MeshedGeometry`
* **`patches`**: A hash map of surface patches, keyed by their topological tag.
* **`edges`**: A hash map of bounding edges/curves, keyed by their topological tag.
* **`edge2faces` & `edge_map`**: Topological mapping structures. `edge_map` maps the tags of a working simulation mesh back to the original baseline tags of this geometric representation (crucial if the simulation mesh undergoes splitting or re-tagging).

#### Workflow & Usage
1. **Initialization**: Instantiated via `MeshedGeometry::new(&boundary_mesh)`. The provided boundary mesh must have correctly tagged faces and edges (e.g., using `mesh.fix()`).
2. **Topology Mapping**: Use `set_topo_map(&mesh_topo)` to reconcile edge tags if the simulation mesh's topology differs slightly from the baseline geometry.
3. **Projection**: When `project()` is invoked, the geometry inspects the dimension of the `TopoTag` to route the point:
   * **`Dim == D`**: Projects the point onto the corresponding 2D surface patch.
   * **`Dim == D - 1`**: Projects the point onto the corresponding 1D boundary curve/edge.
   * **`Dim == 0`**: Represents a hard corner/vertex; the point remains fixed in space.

---

## Utility Functions

* **`curvature(&self, pt: &Vertex<D>, tag: Tag)`**
  Extracts the principal curvature directions at a specific vertex on a surface patch.
* **`write_curvature(&self, fname: &str)`**
  Exports the computed curvature vector fields as `.vtu` files (appending the topological tag to the filename). This is extremely useful for verifying the curvature field visually in tools like ParaView.