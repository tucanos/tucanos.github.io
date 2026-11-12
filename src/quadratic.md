# Quadratic simplices

`tmesh` supports second-order simplex elements:

- quadratic edge (3 nodes),
- quadratic triangle (6 nodes),
- quadratic tetrahedron (10 nodes).

## High-order node layout (quadratic edge/triangle/tetrahedron)

For quadratic simplices in `tmesh`, node sets are split into:

- corner nodes (the linear simplex vertices),
- edge nodes (one node per topological edge),
- no face-interior or volume-interior nodes at order 2.

### Quadratic edge (3 nodes)

Local node ordering is:

- node 0: endpoint $v_0$,
- node 1: endpoint $v_1$,
- node 2: edge node on $(v_0,v_1)$.

So a quadratic edge has 2 corner nodes + 1 edge node.

### Quadratic triangle (6 nodes)

Local node ordering is:

- corner nodes: 0, 1, 2,
- edge nodes:
	- node 3 on edge $(0,1)$,
	- node 4 on edge $(1,2)$,
	- node 5 on edge $(2,0)$.

So a quadratic triangle has 3 corner nodes + 3 edge nodes.

### Quadratic tetrahedron (10 nodes)

Local node ordering is:

- corner nodes: 0, 1, 2, 3,
- edge nodes:
	- node 4 on edge $(0,1)$,
	- node 5 on edge $(1,2)$,
	- node 6 on edge $(2,0)$,
	- node 7 on edge $(0,3)$,
	- node 8 on edge $(1,3)$,
	- node 9 on edge $(2,3)$.

So a quadratic tetrahedron has 4 corner nodes + 6 edge nodes.

### Consistency across neighboring elements

Edge nodes are global mesh vertices. Two neighboring high-order elements that share a linear edge also share exactly the same quadratic edge node index. This is the key continuity assumption used by interpolation, geometry evaluation, and validity checks.

## Linear to quadratic
High-order meshes are generated from linear meshes by inserting one midpoint per unique linear edge.

### Index construction from global edge maps

The construction is deterministic:

1. Build a global edge map from sorted endpoint pairs (`edges()`).
2. Append one new vertex per unique edge (midpoint coordinates).
3. Rebuild each element with corner nodes first, then edge-midpoint nodes using the global edge id map.

For quadratic tetrahedra, local edge-center indices are:

- $(0,1)\rightarrow 4$,
- $(1,2)\rightarrow 5$,
- $(2,0)\rightarrow 6$,
- $(0,3)\rightarrow 7$,
- $(1,3)\rightarrow 8$,
- $(2,3)\rightarrow 9$.

This guarantees that shared edges reuse the same midpoint vertex id.

Rust example:

```rust
use tmesh::mesh::{Mesh, Mesh2d, QuadraticMesh2d, rectangle_mesh, to_quadratic::to_quadratic_triangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 8, 1.0, 8);
let qmsh: QuadraticMesh2d = to_quadratic_triangle_mesh(&msh);
assert!(qmsh.n_verts() > msh.n_verts());
```

Python example:

```python
from pytucanos.mesh import get_square
from pytucanos import Mesh2d

coords, elems, etags, faces, ftags = get_square()
msh = Mesh2d(coords, elems, etags, faces, ftags)
qmsh = msh.to_quadratic()
```

## Quadratic to linear

The inverse conversion is conceptually straightforward: keep only corner-node connectivity and discard edge-node connectivity.

In current APIs, linear-to-quadratic has dedicated helpers, while quadratic-to-linear is typically done explicitly by projecting each quadratic element to its corner indices.

For order-2 simplices, corner-node extraction is:

- quadratic edge $(n_0,n_1,n_2)$ $\rightarrow$ linear edge $(n_0,n_1)$,
- quadratic triangle $(n_0,\dots,n_5)$ $\rightarrow$ linear triangle $(n_0,n_1,n_2)$,
- quadratic tetrahedron $(n_0,\dots,n_9)$ $\rightarrow$ linear tetrahedron $(n_0,n_1,n_2,n_3)$.

## Lagrange versus Bezier control points

Two representations are used in algorithms:

- Lagrange form for mesh storage/interpolation nodes,
- Bezier form for robust geometric bounds.

For a quadratic edge with endpoints $P_0,P_1$ and Lagrange midpoint $M$, the Bezier edge control point is:

$$
B = \frac{1}{2}\left(4M - P_0 - P_1\right).
$$

Equivalent edge-wise formulas are used in triangles and tetrahedra.

In `tmesh`, mesh connectivity and storage are Lagrange-based; Bezier conversion is performed on demand for bound/validity computations.

## Projection and barycentric coordinates

Projection strategy follows simplex recursion:

- if barycentric coordinates are all non-negative, the projection is inside the element,
- otherwise the projection is delegated to boundary faces and the closest face projection is retained.

For quadratic simplices, inverse mapping is nonlinear:

- quadratic edge: stationary distance points from cubic roots,
- quadratic triangle: barycentric coordinates from minimization (Newton-CG),
- quadratic tetrahedron: direct `bcoords` inversion is currently not implemented.

## Validity and distortion via adaptive Jacobian-determinant bounds

For quadratic elements, `tmesh` can bound $\det J$ using Bernstein coefficients:

1. Sample determinant values at Lagrange points.
2. Convert to Bezier coefficients.
3. Recursively subdivide to refine min/max bounds.
4. Detect invalid elements if the lower bound crosses zero.

Using Bernstein basis functions $B_i$:

$$
\det J(\boldsymbol\lambda) = \sum_{i=1}^{N_b} b_i B_i(\boldsymbol\lambda),
$$

with convex-hull bounds:

$$
\min_i b_i \le \det J(\boldsymbol\lambda) \le \max_i b_i.
$$

For quadratic triangles (degree 2, $N_b=6$):

$$
\det J(\boldsymbol\lambda) =
b_1\lambda_1^2 + b_2\lambda_2^2 + b_3\lambda_3^2
+ 2b_4\lambda_1\lambda_2 + 2b_5\lambda_2\lambda_3 + 2b_6\lambda_3\lambda_1.
$$

For quadratic tetrahedra (degree 3, $N_b=20$):

$$
\det J(\boldsymbol\lambda) = \sum_{a_1+a_2+a_3+a_4=3}
\binom{3}{a_1,a_2,a_3,a_4}
b_{a_1a_2a_3a_4}
\lambda_1^{a_1}\lambda_2^{a_2}\lambda_3^{a_3}\lambda_4^{a_4}.
$$

Distortion estimate:

$$
\kappa_J = \frac{d_{\max}}{d_{\min}},
\qquad
d_{\min}=\min_i b_i,
\quad
d_{\max}=\max_i b_i.
$$

Reference:
J.-F. Remacle, J. Lambrechts, B. Seny, E. Marchandise, A. Johnen, C. Geuzaine,
Efficient Visualization of High-Order Finite Elements.
https://gmsh.info/doc/preprints/gmsh_curved_preprint.pdf
