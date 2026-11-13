# Mesh quality

`tmesh` exposes a small set of geometric quality indicators. These indicators are
descriptive: none of them modifies the mesh and none of them returns an error. The only
failure mode is the quadratic-element `gamma` panic listed at the end of this page.

## Face skewness

`face_skewnesses(&all_faces)` iterates over the **internal** faces of the mesh (faces
shared by exactly two elements) and returns `(i0, i1, s)`, where `i0` and `i1` are the
two adjacent elements and $s$ is the skewness of their shared face:

$$
s = \frac{d\left(c_f,\ \mathrm{line}(c_0, c_1)\right)}{|c_1 - c_0|},
$$

with $c_0$ and $c_1$ the centers of the two elements and $c_f$ the center of the shared
face. The denominator is the distance between the two element centers, so $s$ is
dimensionless:

- $s = 0$ when the face center lies on the line joining the element centers, as in an
  orthogonal structured mesh,
- $s$ grows when the face is off-center with respect to the two cells,
- $s$ can exceed $1$ for strongly distorted meshes.

In 2D, a face is an edge and $c_f$ is its midpoint. The `all_faces` argument must be
the face connectivity of the mesh, typically obtained with `Mesh::all_faces()`. Faces
with a number of adjacent elements different from two are skipped.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 3, 1.0, 3);
let all_faces = msh.all_faces();

let skewnesses: Vec<_> = msh.face_skewnesses(&all_faces).collect();

assert_eq!(skewnesses.len(), 8); // internal edges of the 2 x 2 cell grid
assert!(skewnesses.iter().all(|&(_, _, s)| s.abs() < 1e-6));
```

## Edge length ratio

`edge_length_ratios()` returns, for every element, the ratio between its longest and
its shortest edge:

$$
r_K = \frac{\max_{e \subset K} |e|}{\min_{e \subset K} |e|} \geq 1.
$$

The ratio is $1$ for an equilateral simplex and $\sqrt{2}$ for a right isosceles
triangle. It detects stretched elements, but not all bad elements: a flat tetrahedron
can have almost equal edge lengths.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 3, 1.0, 3);

// Every triangle of the structured grid is right isosceles.
let sqrt2 = 2.0_f64.sqrt();
assert!(msh
	.edge_length_ratios()
	.all(|r| (r - sqrt2).abs() < 1e-12));
```

## Element quality

`elem_gammas()` returns the geometric quality $\gamma$ of every element, normalized so
that the regular simplex has $\gamma = 1$ and degenerate elements have
$\gamma \to 0$:

| Element | $\gamma$ |
| --- | --- |
| `Edge` | $1$ (always) |
| `Triangle` | $\dfrac{4 \sin A \sin B \sin C}{\sin A + \sin B + \sin C} = \dfrac{2 r_{\mathrm{in}}}{R}$ |
| `Tetrahedron` | $\dfrac{3 r_{\mathrm{in}}}{R}$ |

Here $A, B, C$ are the triangle angles, $r_{\mathrm{in}}$ is the inscribed radius, and
$R$ is the circumradius. For a tetrahedron, $r_{\mathrm{in}} = 3V/S$ where $V$ is the
volume and $S$ the total area of the four faces. The circumradius computation is
ill-conditioned for nearly coplanar vertices; the implementation returns $0$ in that
case (and for zero-volume tetrahedra).

This indicator is useful to locate slivers and needles, especially in 3D, where the
edge ratio alone does not see flat elements. It is only implemented for linear
elements: quadratic triangles and tetrahedra do not provide a `gamma` and trigger an
`unreachable` panic.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 3, 1.0, 3);

// gamma = 2 (sqrt(2) - 1) for a right isosceles triangle.
let expected = 2.0 * (2.0_f64.sqrt() - 1.0);
assert!(msh
	.elem_gammas()
	.all(|g| (g - expected).abs() < 1e-12));
```

Integration and $L^2$ norms are documented separately in
[Integration](integration.md).

> [!TIP]
> **Quality panics**
>
> | Message | Meaning / action |
> | --- | --- |
> | `internal error: entered unreachable code` | `elem_gammas` was called on a quadratic mesh; `gamma` is only implemented for linear triangles and tetrahedra. |
>
> `face_skewnesses`, `edge_length_ratios`, and `elem_gammas` return iterators and do not
> allocate a result array; they are usually consumed with `collect`, `map`, `min`, or
> `max`. The face skewness requires the full face connectivity of the mesh, not a
> filtered or submesh connectivity.
