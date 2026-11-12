# Extrusion

Extrusion turns a 2D triangle mesh into a single layer of 3D prisms. It is a cheap way
to obtain a 3D test case, to build a prismatic mesh from a planar one, or to visualize
a 2D field as a thin 3D volume.

`ExtrudedMesh2d` is not a `Mesh<3>`: it stores prisms, triangular caps, and
quadrilateral sides separately, and provides its own accessors. The only file format it
can write is `.vtu` (`write_vtk`); use `to_mesh2d` to go back to a 2D mesh (see below).

## Construction

- `ExtrudedMesh2d::from_mesh2d(&msh, h)` extrudes a `Mesh2d` by a signed distance `h`.
- `msh.extrude::<T>(h)` is the same operation as a method on the 2D mesh.
- `ExtrudedMesh2d::new(...)` builds one explicitly from coordinates, prisms, triangles,
  quads, and their tags.

Let $n$ be the number of vertices of the 2D mesh. The extruded mesh contains $2n$
vertices: the layer at $z = 0$ keeps the original numbering $0 \ldots n-1$, and the
layer at $z = h$ uses $i + n$. Every triangle $(v_0, v_1, v_2)$ becomes one prism
joining the two layers.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 3, 2.0, 4);

let extruded = msh.extrude::<usize>(0.5);

assert_eq!(extruded.n_verts(), 2 * msh.n_verts()); // 24
assert_eq!(extruded.n_prisms(), msh.n_elems()); // 12
assert_eq!(extruded.n_tris(), 2 * msh.n_elems()); // 24 caps
assert_eq!(extruded.n_quads(), msh.n_faces()); // 10 sides
```

`pytucanos` exposes the same operation:

```python
import numpy as np
from pytucanos import Mesh2d

msh = Mesh2d.rectangle_mesh(np.linspace(0.0, 1.0, 3), np.linspace(0.0, 2.0, 4))
extruded = msh.extrude(0.5)

assert extruded.n_verts() == 24
assert extruded.n_prisms() == 12
assert extruded.n_tris() == 24
assert extruded.n_quads() == 10
```

## Orientation

The prism local numbering depends on the sign of `h`:

- for $h > 0$, the prism is $(v_0, v_1, v_2, v_0+n, v_1+n, v_2+n)$,
- for $h < 0$, the second and third vertices of each layer are swapped,
  $(v_0, v_2, v_1, v_0+n, v_2+n, v_1+n)$.

This keeps the signed volume of every prism positive regardless of the extrusion
direction. The top cap is ordered with the opposite orientation of the bottom cap, so
the two ends are consistently oriented.

The construction does not validate that the 2D triangles are positively oriented. If in
doubt, fix the 2D mesh first (see [Validity](validity.md)).

## Tags

Extrusion creates new boundary surfaces that have no counterpart in the 2D mesh, so
they receive dedicated tags:

| Entity | Tag |
| --- | --- |
| Prisms | inherited from the 2D element tags `etags` |
| Side quadrangles | inherited from the 2D face tags `ftags` |
| Bottom cap triangles ($z = 0$) | `Tag::MAX` |
| Top cap triangles ($z = h$) | `Tag::MAX - 1` |

With the default `i16` tag type, `Tag::MAX = 32767` and `Tag::MAX - 1 = 32766`. These
values change with the `32bit-tags` and `64bit-tags` features (see [Tags](tags.md)).

Negative 2D face tags, which mark protected boundaries in remeshing workflows, are
preserved on the corresponding side quadrangles.

## Unextrusion

`ExtrudedMesh2d::to_mesh2d()` extracts the $z = 0$ layer and rebuilds a `Mesh2d`:

- the first $n$ vertices are projected to the plane,
- each prism $(v_0, v_1, v_2, \ldots)$ becomes the triangle $(v_0, v_1, v_2)$,
- each side quadrangle becomes the edge joining its first two vertices,
- element and face tags are copied back,
- element and face orientations are fixed.

The result is independent of the sign of `h`: it always represents the bottom layer.
For $h > 0$, extruding and unextruding returns the original mesh.

> [!WARNING]
> For $h < 0$, the reconstructed element connectivity is a cyclic rotation of the input
> one, so the mesh is geometrically identical but its arrays differ. `check_equals`
> reports a mismatch such as `Triangle([0, 1, 6]) != Triangle([6, 0, 1])` even though
> both triangles are the same. Compare volumes or reorder the vertices if an exact array
> match is required.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 10, 2.0, 15);
let extruded = msh.extrude::<usize>(1.0);

let back: Mesh2d = extruded.to_mesh2d();
msh.check_equals(&back, 1e-12).unwrap();
```

## Extruding a dual mesh

`DualMesh2d::extrude(h)` extrudes a vertex-centered dual mesh into a
`SimplePolyMesh<3>`: one polyhedral cell per 2D dual polygon, with prismatic side
faces. It is used to build 3D control volumes from a 2D dual mesh.

Each 2D dual face must be an edge (two vertices), and the boundary of each dual element
must merge into a single closed polygon. The caps use the same `Tag::MAX` (bottom) and
`Tag::MAX - 1` (top) convention, and the lateral faces keep the 2D dual face tags.

```rust
# extern crate tmesh;
use tmesh::dual::{DualMesh, DualMesh2d, DualType, PolyMesh};
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 10, 2.0, 15);
let dual = DualMesh2d::new(&msh, DualType::Median);

let extruded = dual.extrude(1.0);
assert_eq!(extruded.n_elems(), dual.n_elems());
```

## Errors and panics

> [!TIP]
> **Extrusion panic messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `assertion failed: v[2].abs() < 1e-12` | `to_mesh2d` found a vertex of the bottom layer that is not in the plane $z = 0$. The extruded mesh was probably built manually with a shifted first layer. |
> | ``assertion `left == right` failed`` with `right: 2` | `DualMesh2d::extrude` requires every 2D dual face to be an edge. The dual mesh was built with non-edge faces. |
> | ``assertion `left == right` failed`` with `right: 1` | The boundary of a 2D dual element did not merge into a single closed polygon. The dual connectivity is inconsistent. |
>
> No other validation is performed: `h = 0` is accepted and produces degenerate
> prisms, and no check is made on the orientation of the input triangles.
