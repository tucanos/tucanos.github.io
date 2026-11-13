# Dual meshes

This page documents the assumptions and algorithms used to build vertex-centered dual meshes in `tucanos` (implemented in `tmesh`).

## Goal

Given a primal simplex mesh:

- triangle mesh in 2D,
- tetrahedral mesh in 3D,

build a dual control-volume structure centered on primal vertices.

`tmesh` computes two consistent outputs:

- an edge-based representation (edge-associated normals and boundary-face normals),
- an explicit polygonal/polyhedral dual mesh.

## Inputs and outputs

Inputs:

- primal simplex mesh (geometry + connectivity),
- face tags (for boundary contributions),
- dual center policy (`Median`, `Barth`, or `ThresholdBarth(\ell)`).

Outputs:

- dual vertices,
- dual faces,
- dual elements (one dual element per primal vertex),
- edge-associated dual normals,
- boundary dual-face normals/tags.

## Dual center policies

`tmesh` supports three center policies:

- `Median`: geometric centers of simplices,
- `Barth`: circumcenter-based centers,
- `ThresholdBarth(\ell)`: Barth center with barycentric clipping threshold $\ell$.

In obtuse or near-degenerate cases, Barth-like constructions can fall back to lower-dimensional centers to preserve robust dual-cell construction.

## 2D algorithm (triangular primal mesh)

For a primal triangle mesh, the dual construction proceeds as follows:

1. Enumerate primal edges and boundary vertices.
2. Create dual vertices from:
	- boundary primal vertices,
	- primal-edge midpoints,
	- one center per primal triangle (depends on policy).
3. Build internal dual segments linking edge midpoints to triangle centers.
4. Build boundary dual segments linking boundary vertices to boundary-edge midpoints.
5. Assemble oriented face incidences around each primal vertex to define dual polygons.
6. Accumulate oriented edge normals from dual segments.
7. Remove empty/degenerate dual faces and zero-length normal contributions.

The result is a polygonal dual mesh with one dual cell per primal vertex.

## 3D algorithm (tetrahedral primal mesh)

For a primal tetrahedral mesh, the dual construction proceeds as follows:

1. Enumerate primal edges and primal faces.
2. Create dual vertices from:
	- boundary primal vertices,
	- edge midpoints,
	- face centers,
	- tetrahedron centers.
3. Build internal dual triangles from edge-center/face-center/cell-center triplets.
4. Build boundary dual triangles from boundary-vertex/edge-center/face-center triplets.
5. Deduplicate dual faces (required for Barth-like center policies).
6. Assemble oriented face incidences per primal vertex to define dual polyhedra.
7. Accumulate and filter edge-based dual normals.

Again, one dual cell is created per primal vertex.

## Invariants and interpretation

Important structural properties:

- dual element count equals primal vertex count,
- dual volumes partition primal volume (up to floating-point tolerance),
- boundary dual faces inherit boundary tags from primal boundary topology,
- edge normals are oriented and scaled for FV-like flux integrations.

## Validity checks and errors

`PolyMesh::check()` (and the dimension-specific helper behind it) validates the dual
connectivity: consistent array lengths, in-range indices, closed elements, and unique
faces referenced at most once per element with opposite orientations.

> [!TIP]
> **`PolyMesh::check()` error messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `Inconsistent sizes (faces)` | The number of face tags differs from the number of faces. |
> | `Inconsistent indices (faces)` | A dual face references a vertex index outside the dual vertex array. |
> | `Inconsistent indices (elems)` | A dual element references a face index outside the dual face array. |
> | `Element {i} invalid: vol={v} < 0 ({e:?})` | A dual element has non-positive volume. The dual construction is degenerate, usually because the primal mesh has a flat or inverted element. |
> | `Element {i} not closed ({e:?})` | The face incidences of a dual element do not form a closed boundary. |
> | `Face {i} ({f:?}) = face {j} ({f_j:?})` | The same dual face appears twice in the face list. |
> | `Face {i} appears multiple times in element {i_elem}` | A dual element lists the same face more than once. |
> | `Face {i} appears twice with True` / `Face {i} appears twice with False` | A face is referenced by two elements with the same orientation. Orient one of them oppositely. |
> | `Face {i} appears 3 times` | A dual face is shared by more than two elements (non-manifold dual mesh). |
> | `Check only available with simplex faces` | `check()` supports only `PolyFaceType::Simplices`; polygonal/polyhedral face storage has no generic check. |
>
> Polygon assembly can also fail during construction with `unable to merge faces`
> when a set of candidate faces cannot be merged into a single face.

> [!TIP]
> **Dual-mesh panics**
>
> Several dual-mesh accessors require simplex faces and **panic** instead of returning
> an error: `Faces are not simplices` (volume, first moment, and tagged-normal
> computations), `Element {i} has zero volume` (first-moment computation),
> `No free slot for dual face insertion` / `No free slot for internal
> dual face insertion` / `No free slot for boundary dual face insertion` (internal
> insertion invariant), and `is_vertex_in_elem is only implemented for D=2 and D=3`.
> These indicate a misuse of the API or an internal invariant violation, not invalid
> input data.

## Rust example

```rust
# extern crate tmesh;
use tmesh::{
	 dual::{DualMesh, DualMesh2d, DualType, PolyMesh},
	 mesh::{Mesh, Mesh2d, rectangle_mesh},
};

let msh = rectangle_mesh::<Mesh2d>(1.0, 10, 2.0, 15);
let dual = DualMesh2d::new(&msh, DualType::Median);
dual.check().unwrap();

assert_eq!(dual.n_elems(), msh.n_verts());
```

## Python example (`pytucanos`)

```python
import numpy as np
from pytucanos import Mesh2d, DualMesh2d, DualType

x = np.linspace(0.0, 1.0, 10)
y = np.linspace(0.0, 2.0, 15)
msh = Mesh2d.rectangle_mesh(x, y)

dual = DualMesh2d(msh, DualType.Median)
assert dual.n_elems() == msh.n_verts()

dual_bdy, _ = dual.boundary()
dual_bdy.fix()
```
