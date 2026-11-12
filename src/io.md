# Mesh I/O

Mesh input/output is provided by the `Mesh` trait (`tmesh::mesh`) and by the
`tmesh::io` module. The formats used in practice are:

- `.meshb` for meshes and `.solb` for fields (through `minimeshb`),
- `.vtu` (VTK XML unstructured grid) for visualization and post-processing.

## `.meshb` meshes

`Mesh::from_meshb(file_name)` reads a `.meshb` mesh and `Mesh::write_meshb(file_name)`
writes one. Linear and quadratic edges, triangles, and tetrahedra are supported.
Missing sections are allowed when reading (for example a mesh without a face section).
Element and face tags are preserved.

The tag width is selected at compile time: `Tag` is `i16` by default, `i32` with the
`32bit-tags` feature, and `i64` with the `64bit-tags` feature.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 10, 2.0, 15);
msh.write_meshb("mesh.meshb").unwrap();

let msh2 = Mesh2d::from_meshb("mesh.meshb").unwrap();
assert_eq!(msh2.n_verts(), msh.n_verts());
```

> [!TIP]
> **`.meshb` error and warning messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `Tag {t} is too large for the tag type` | A tag read from a `.meshb` file does not fit the configured `Tag` type (`i16` by default). Reduce the tag value or enable the `32bit-tags` / `64bit-tags` feature. |
> | `not reading faces when elements are edges` | Warning: the element type is a 1D edge, so the face section is ignored. |
> | `skip faces in meshb export` | Warning: the mesh stores 0D faces (nodes), which `.meshb` output does not write. |
> | `not yet implemented` | Panic: the simplex order (only 1 and 2) or the number of vertices is not supported by the `.meshb` reader/writer. |
>
> With `64bit-tags`, writing converts each tag to the 32-bit integer stored by `.meshb`
> and panics if a tag does not fit in `i32`. Errors reported by `minimeshb` itself
> (invalid file, unsupported version, truncated sections, I/O failures) are propagated
> unchanged.

## `.solb` solution fields

`Mesh::write_solb(values, file_name, location)` writes a field and
`Mesh::read_solb(file_name)` reads one. `location` selects where the values live:

- `SolutionLocation::Vertices`,
- `SolutionLocation::Elements`,
- `SolutionLocation::Faces`,
- `SolutionLocation::Edges`.

Values are stored component-interleaved. For vertex, element, and face locations the
array length must be `n_components * n_entities`; edge locations are treated as scalar
data (see the warning below). Reading returns the flattened values and the number of
components per entity; only vertex solutions (`SolAtVertices`) are read.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, SolutionLocation, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 10, 2.0, 15);
let field: Vec<f64> = msh.verts().map(|x| x[0] + x[1]).collect();

msh.write_solb(&field, "field.solb", SolutionLocation::Vertices).unwrap();
let (values, n_comp) = Mesh2d::read_solb("field.solb").unwrap();
assert_eq!(n_comp, 1);
assert_eq!(values.len(), msh.n_verts());
```

> [!TIP]
> **`.solb` error and panic messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `Invalid dimension {d} in {file_name} (expected {D})` | `read_solb` opened a solution whose spatial dimension differs from the mesh dimension. |
> | `Invalid solution size {len} for {n} entities` | The array passed to `write_solb` is empty or its length is not a multiple of the number of target entities. |
> | `Solutions with {n_comp} components are not supported in {D}D` | `write_solb` accepts 1, 2, or 3 components in 2D and 1, 3, or 6 components in 3D. |
> | `Solutions with {m} components are not supported in {d}D` | `read_solb` found an unsupported component count in the file (same allowed values as above). |
> | `Solution at elements is only implemented for linear simplicial meshes` | Panic: element-located solutions require linear triangles or tetrahedra. |
> | `Solution at faces is only implemented for linear simplicial meshes` | Panic: face-located solutions require linear edges, triangles, or tetrahedra. |
> | `Solution at edges is only implemented for linear simplicial meshes` | Panic: edge-located solutions require a linear mesh. |

> [!WARNING]
> `write_solb` with `SolutionLocation::Edges` assumes scalar data and does not verify
> that the array length matches the number of edges. `read_solb` only reads vertex
> solutions; `.solb` files storing element, face, or edge solutions cannot be read back
> with it.

## `.vtu` output

`Mesh::write_vtk(file_name)` writes an XML `.vtu` file and returns a
`std::io::Result`. Element tags are written as the cell data array `tags`. Points of
2D meshes are lifted to `z = 0`.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 10, 2.0, 15);
msh.write_vtk("mesh.vtu").unwrap();
```

For more control, use `tmesh::io::VTUFile` directly:

```rust
# extern crate tmesh;
use tmesh::io::VTUFile;
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 10, 2.0, 15);

let mut writer = VTUFile::from_mesh(&msh);
writer.add_point_data("x", 1, msh.verts().map(|v| v[0]));
writer.add_cell_data("etag", 1, msh.etags().map(|t| t as f64));
writer.export("mesh.vtu").unwrap();
```

`VTUFile` also provides `from_extruded_mesh(...)` for `ExtrudedMesh2d` and
`from_poly_mesh(...)` for dual `PolyMesh` structures. Cell types are written as:

| Mesh elements | VTK cell type |
| --- | --- |
| linear edges, triangles, tetrahedra | 3, 5, 10 |
| quadratic edges, triangles, tetrahedra | 21, 22, 24 |
| prisms (`ExtrudedMesh2d`) | 13 |
| polygons, polyhedra (`PolyMesh`) | 7, 42 |

> [!TIP]
> **`.vtu` error and panic messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `faces = {faces:?}` | Internal assertion while assembling a polygonal `PolyMesh`: the face loop did not merge into a single polygon. |
> | `not yet implemented` | `PolyMeshType::Polylines` export is not implemented; convert the polyline mesh to polygons or write it another way. The same panic occurs for meshes with dimension other than 2 or 3, and for unsupported simplex orders/types. |
>
> `export` returns the underlying `std::io::Result` for file creation and write
> failures; the VTK writer itself does not return crate-level validation errors.

## STL input

Surface meshes can also be read from STL files with `tmesh::mesh::read_stl`. The
reader merges vertices with identical bit patterns and tags all triangles with `1`; it
does not create boundary faces or fix the facet orientation. See
[Mesh generation](mesh_generation.md#reading-stl-surfaces) for the details and for an
example.
