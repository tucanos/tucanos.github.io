# Element types

`tmesh` stores elements as topological simplices. The mesh type is parameterized by an ambient dimension $D$ and an element type `C` implementing the `Simplex` trait, e.g. `GenericMesh<2, Triangle<usize>>`.

## Linear simplices

| Element | Dim. | Nodes | Typical use |
| --- | --- | --- | --- |
| `Edge` | 1 | 2 | 1D mesh, boundary/interface edges |
| `Triangle` | 2 | 3 | 2D mesh, boundary/interface triangles |
| `Tetrahedron` | 3 | 4 | 3D volume mesh |

The faces of an element are also simplices: edges for a triangle, triangles for a tetrahedron, and nodes for an edge.

Common aliases:

- `Mesh2d` = `GenericMesh<2, Triangle<usize>>`,
- `Mesh3d` = `GenericMesh<3, Tetrahedron<usize>>`,
- `BoundaryMesh2d` = `GenericMesh<2, Edge<usize>>`,
- `BoundaryMesh3d` = `GenericMesh<3, Triangle<usize>>`.

## Quadratic simplices

Second-order elements have corner nodes plus one node per topological edge:

| Element | Nodes |
| --- | --- |
| `QuadraticEdge` | 3 |
| `QuadraticTriangle` | 6 |
| `QuadraticTetrahedron` | 10 |

See [Quadratic simplices](quadratic.md) for the local node ordering.

Aliases mirror the linear case: `QuadraticMesh2d`, `QuadraticMesh3d`, `QuadraticBoundaryMesh2d`, and `QuadraticBoundaryMesh3d`.

## Canonical volume cells

Non-simplex cells can be inserted into a simplex mesh; they are subdivided deterministically into simplices:

| Cell | Nodes | Subdivision |
| --- | --- | --- |
| `Quadrangle` | 4 | 2 triangles |
| `Pyramid` | 5 | 2 tetrahedra |
| `Prism` | 6 | 3 tetrahedra |
| `Hexahedron` | 8 | 5 or 6 tetrahedra |

See [Subdivision of canonical elements](subdivision.md) for the exact rules and the limitation to linear cells.

## Working with elements

- [Simplices](simplices.md): topological versus geometric simplices, local numbering, and geometric operations.
- [Tags](tags.md): element and face tags.
- [Adjacencies](adjacencies.md): connectivities derived from element and face arrays.
