# Isosurface splitting

`Mesh::split_isosurface(f)` splits a mesh along the `0.0` level set of a scalar field
`f` defined at the mesh vertices. The two sides of the level set are no longer joined
by volume elements; they meet only at the new explicit interface. This is useful to
insert an isosurface into a volume mesh, to separate regions before adaptation, or to
extract an interface.

The input must be a linear simplex mesh (`order() == 1`) whose element tags are all `1`
and whose face tags are all positive.

## Output conventions

After splitting:

- elements with `f > 0` are tagged `1`,
- elements with `f < 0` are tagged `-1`,
- existing face tags are kept on the positive side and negated on the negative side,
- the new interface faces are tagged `Tag::MAX`.

`Tag::MAX` is `i16::MAX` (32767) by default, or the maximum of the configured tag type.
Because negative face tags are treated as protected by the remesher, the flipped tags
on the negative side preserve the original boundary protections.

The method returns `(new_mesh, split_edges)`, where `split_edges` maps every original
edge cut by the isosurface to the index of the vertex inserted on that edge.

## Algorithm

1. Validate the field and replace exact zeros by `1e-12 * max|f|` to avoid
   degenerate intersections.
2. For every edge with `f[i] * f[j] <= 0`, insert a vertex at
   $t = -f_i / (f_j - f_i)$ along the edge.
3. Split each element according to the signs of its vertex values:
   - **1D edge**: split into two edges.
   - **2D triangle**: one triangle on the single-vertex side and a quad on the other
    side; the quad is split into two triangles using the global-index rule.
   - **3D tetrahedron**: a `1 + 3` split produces one tetrahedron and one prism
     (decomposed into three tetrahedra), while a `2 + 2` split produces one prism per
     side. The interface is a triangle for `1 + 3` and a quad (two triangles) for
     `2 + 2`.
4. Split the existing boundary/interface faces the same way, flipping their tags on
   the negative side. Faces of 2D meshes are edges and faces of 3D meshes are
   triangles; 0D faces (nodes) are not supported.
5. Add the new interface faces with tag `Tag::MAX`.

The decomposition of prisms and quads uses the same global-index rule as the
subdivision routines, so the result is independent of the local element ordering.

## Example

```rust
# extern crate tmesh;
use tmesh::mesh::{BoundaryMesh2d, Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 20, 1.0, 20);
let f: Vec<f64> = msh.verts().map(|x| x[0] - 0.5).collect();

let (split, _split_edges): (Mesh2d, _) = msh.split_isosurface(&f).unwrap();

// Positive and negative regions.
assert!(split.etags().all(|t| t == 1 || t == -1));

// Extract the interface (Tag::MAX) as a boundary mesh of edges.
let (interface, _): (BoundaryMesh2d, _) =
	split.extract_faces(|t| t == tmesh::Tag::MAX);
assert!(interface.n_elems() > 0);
```

## Errors and panics

> [!TIP]
> **`split_isosurface()` error and panic messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `Invalid field size {len} (expected {n})` | `f` must provide exactly one value per vertex. |
> | `The field contains non finite values` | `f` contains `NaN` or `Inf`; clean the field first. |
> | `The field is identically zero` | The level set is the whole mesh; provide a non-trivial field. |
> | `Exact zeros encountered. Ensure f is perturbed slightly before the element loop.` | Internal panic: exact zeros were not removed before element splitting. This should not occur through the public API; report it as a bug. |
> | `not yet implemented` | Panic: an unsupported element or face type was reached. Elements must be linear edges, triangles, or tetrahedra; faces must be linear edges or triangles. |
>
> The function also asserts that the element tags are all `1`, that all face tags are
> positive, and that elements are linear (`order() == 1`).
