# Validity

A mesh is considered valid if
- all vertices are used by at least one element,
- all listed faces belong to at least one element,
- all required tagged faces exist:
	- external boundary faces,
	- interfaces between different element tags,
	- non-manifold interfaces,
- tagged internal faces are not incorrectly placed inside a same-tag region.
- all element signed measures are non-negative,
 

## What `check()` does

The `check` methods checks the validity criteria mentioned above as well as that
- sizes are consistent:
	- number of element tags equals number of elements,
	- number of face tags equals number of faces,
- all element/face indices are in bounds,
- boundary face orientation is valid when orientation is defined for the mesh type,
- global volume consistency checks pass.

> [!TIP]
> **`Mesh::check()` error messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `Inconsistent sizes (elems)` | The number of element tags differs from the number of elements. |
> | `Inconsistent sizes (faces)` | The number of face tags differs from the number of faces. |
> | `{n} vertices are not used in any element` | Unused vertices are present; remove them or fix the connectivity. |
> | `Invalid index in elems` | An element references a vertex index outside the vertex array. |
> | `Elem has a <0 volume` | An element has a negative signed measure. Reorient it, or call `fix()`. |
> | `Invalid index in faces` | A face references a vertex index outside the vertex array. |
> | `Face belong to no element` | A listed face is not a face of any element. |
> | `Invalid face orientation: center = ..., normal = ..., face = ..., tag = ...` | A tagged boundary face is oriented inward. `fix()` reorients boundary faces. |
> | `Invalid volume : {vol} from elements, {vol2} from boundary faces` | The volume computed from elements differs from the volume computed from boundary faces (linear, oriented meshes). Check boundary-face orientation and completeness. |
>
> Missing face tags are reported by `check()` with the messages listed in [Tags](tags.md#face-tagging-errors). `fix()` adds the missing boundary and interface tags before running these checks again.

## Comparing two meshes

`check_equals(other, tol)` verifies that two meshes have the same vertex coordinates
(within `tol`), connectivities, and tags.

> [!TIP]
> **`check_equals()` error messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `Vertex {i}: {v0:?} != {v1:?}` | Vertex `i` differs by more than `tol`. |
> | `Element {i}: {e0:?} != {e1:?}` | Element connectivity `i` differs. |
> | `Element tag {i}: {t0:?} != {t1:?}` | Element tag `i` differs. |
> | `Face {i}: {e0:?} != {e1:?}` | Face connectivity `i` differs. |
> | `Face tag {i}: {t0:?} != {t1:?}` | Face tag `i` differs. |
>
> The comparison is order-sensitive: reordering vertices, elements, or faces must be
> undone (for example with the reordering helpers) before calling `check_equals`.

## What `fix()` does

`fix()` performs an ordered repair pipeline:

1. Reorient elements with negative signed measure.
2. Build full face-to-element connectivity (`all_faces`).
3. Add missing boundary face tags.
4. Add missing internal-interface tags.
5. Reorient boundary/interface faces consistently.
6. Run full validity checks.

Return value:

- generated boundary-tag map,
- generated internal-interface-tag map.

If any required condition still fails after these repair steps, `fix()` returns an error.

## Example: one triangle in 3D with invalid/incomplete boundary faces

The example below builds a `BoundaryMesh3d` with:

- one triangle element,
- an intentionally reversed boundary edge,
- missing boundary edges (only 1 of the 3 is provided).

Calling `fix()` reorients listed boundary faces as needed and adds the missing boundary faces.

```rust
# extern crate tmesh;
use tmesh::{Tag, Vert3d};
use tmesh::mesh::{BoundaryMesh3d, Edge, Mesh, Triangle};

let verts = vec![
	Vert3d::new(0.0, 0.0, 0.0),
	Vert3d::new(1.0, 0.0, 0.0),
	Vert3d::new(0.0, 1.0, 0.0),
];

// One triangle in 3D.
let elems = vec![Triangle::new(0, 1, 2)];
let etags: Vec<Tag> = vec![1];

// Incomplete + intentionally reversed boundary-face list.
// For triangle elements, topological faces are edges.
let faces = vec![Edge::new(1, 0)];
let ftags: Vec<Tag> = vec![10];

let mut msh = BoundaryMesh3d::from_vecs(verts, elems, etags, faces, ftags);

let (_new_bdy_tags, _new_ifc_tags) = msh.fix().unwrap();

// After fix(), the triangle boundary is complete: 3 tagged edges.
assert_eq!(msh.faces().count(), 3);

// Explicit validity check (optional after fix).
let all_faces = msh.all_faces();
msh.check(&all_faces).unwrap();
```

Python (`pytucanos`):

```python
import numpy as np
from pytucanos import BoundaryMesh3d

coords = np.array([
	[0.0, 0.0, 0.0],
	[1.0, 0.0, 0.0],
	[0.0, 1.0, 0.0],
], dtype=float)

# One triangle in 3D.
elems = np.array([[0, 1, 2]], dtype=np.int64)
etags = np.array([1], dtype=np.int16)

# Incomplete + intentionally reversed boundary-face list.
# For triangle elements, topological faces are edges.
faces = np.array([[1, 0]], dtype=np.int64)
ftags = np.array([10], dtype=np.int16)

msh = BoundaryMesh3d(coords, elems, etags, faces, ftags)

bdy_tags, ifc_tags = msh.fix()

# After fix(), the triangle boundary is complete: 3 tagged edges.
assert msh.n_faces() == 3

# Explicit validity check (optional after fix).
msh.check()
```

## Remesher consistency errors

`Remesher::check()` validates the internal remesher data (elements, vertices, edges, tagged faces, and metrics) before and after local operations.

> [!TIP]
> **`Remesher::check()` error messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `Vertex not found` | An element references a vertex that is not stored in the remesher. |
> | `Invalid vertex to element (missing vertex)` | The element is missing from a vertex's element list. |
> | `Invalid vertex to element (missing element)` | A vertex references an element that is not stored in the remesher. |
> | `Negative volume` | An element has negative volume; reorient it. |
> | `Missing edge` | An edge of an element is absent from the remesher edge table. |
> | `Invalid edge (missing vertex)` | An edge references a vertex that is not stored in the remesher. |
> | `Invalid edge {edg:?}: counted in {count} elements, found in {n}` | The edge counter disagrees with the actual number of incident elements. |
> | `Invalid tagged face {face:?} - tag={tag}` | A tagged face references a vertex that is not stored in the remesher. |
> | `A face belongs to more than 2 elements` | The remesher does not accept a non-manifold interior face in the mesh it operates on. Prepare a manifold cavity or tag the interface explicitly. |
> | `A face ... belonging to 2 elements ... with different tags is not tagged correctly` | An internal interface between different element tags is missing its tag. |
> | `A face ... belonging to 2 elements ... with the same tags is not tagged correctly` | An internal face between two same-tag elements is tagged; an internal face inside a same-tag region must not be tagged. |
> | `A face ... belonging to 1 element ... is not tagged correctly` | A boundary face is missing its tag. |
>
> `Remesher::check()` also calls `Metric::check()` on every vertex, so [metric validity errors](metric.md#metric-validity-and-construction) can be returned as well.

> [!TIP]
> **`Remesher::new()`, mutation, and parallel errors**
>
> | Message | Meaning / action |
> | --- | --- |
> | `Vertex not present` | `remove_vertex` was given an index that does not exist. |
> | `Vertex used` | `remove_vertex` was given a vertex still referenced by elements. |
> | `Element not present` | `remove_elem` was given an index that does not exist. |
> | `Element vertex not present: {el:?}` | `insert_elem` references a vertex that has not been inserted. |
> | `Invalid element: {el:?} / {ge:?}, quality: {q}` | `insert_elem` was given an element with non-positive metric quality (flat or inverted). |
> | `Tagged face {face:?} already present` | `Remesher::new` panics when the mesh lists the same tagged face twice. |
> | `At least a vertex of face {face:?} is not in the mesh` | `Remesher::new` panics when a tagged face references an unknown vertex. |
> | `Tagged face not present` | `remove_tagged_face` was given a face that is not in the tagged-face list. |
> | `Failed to remesh partition {i_part}.\n{error}` | `ParallelRemesher` wraps the failure of one partition; the inner message identifies the actual error. |
> | `All the elements are in the interface` | Warning: after remeshing a partition, no face separates the elements tagged `1` from those tagged `2`, so no interface boundary tag is assigned. |
>
> `Remesher::new` also propagates `Geometry::check(topo)` and `insert_elem` errors, and asserts `metric.len() == mesh.n_verts()` and `topo.vtags().len() == mesh.n_verts()`.

