# Graphs and partitioning

`tmesh` uses compressed sparse row (CSR) graphs for connectivities that cannot be stored
as a fixed-width array: vertex-to-vertex, vertex-to-element, and element-to-element
adjacency. The same structure supports reordering and mesh partitioning.

## CSR graphs

A `CSRGraph` stores:

- `ptr`: $n+1$ offsets into the column array,
- `indices`: the column indices, row by row,
- an optional `values` array attached to the column indices,
- `m`: the number of columns.

The neighbours of vertex `i` are `indices[ptr[i]..ptr[i+1]]`. Rows are sorted and
contain no duplicate column. Because the storage is symmetric for adjacency graphs,
each undirected edge appears twice: `n_edges()` is twice the number of edges.

`CSRGraph::new(ptr, indices)` builds a graph from the raw CSR arrays. It is low-level:
the `ptr` array must be consistent, and `indices` must not be empty.

`CSRGraph::from_edges(edges, n_verts)` builds the symmetric graph of an edge iterator,
optionally fixing the number of vertices. Duplicate edges and self-loops are not
allowed: they create duplicate columns and trip an internal assertion.

`CSRGraph::transpose(elems, n_verts)` builds the vertex-to-element graph from an
element-to-vertex connectivity. Its `values` store the local position of the vertex in
the element, which is what `row_and_values` returns. This is the graph used by the
$L^2$-projection gradient and Hessian.

| Method | Meaning |
| --- | --- |
| `n()` | number of rows (vertices) |
| `m()` | number of columns |
| `n_edges()` | number of stored entries |
| `row(i)`, `row_ptr(i)`, `rows()` | neighbours of a vertex |
| `row_and_values(i)` | neighbours and attached values (panics when there are no values) |
| `subgraph(ids)` | induced subgraph, renumbered in the order of `ids` |
| `connected_components()` | component label per vertex |
| `reverse_cuthill_mckee()` | RCM permutation |

`tmesh::graph::reindex(elems)` is a helper that renumbers the vertex indices of an
element array by order of first appearance and returns the old-to-new map.

```rust
# extern crate tmesh;
use tmesh::graph::CSRGraph;

// Two components: the path 0 - 1 - 2 and the edge 3 - 4.
let graph =
	CSRGraph::from_edges([[0, 1], [1, 2], [3, 4]].into_iter(), None);

assert_eq!(graph.n(), 5);
assert_eq!(graph.n_edges(), 6); // two entries per undirected edge
assert_eq!(graph.row(1), [0, 2]);

// The component labels follow the order in which components are discovered.
assert_eq!(graph.connected_components().unwrap(), [0, 0, 0, 1, 1]);

// The reversed Cuthill-McKee order.
assert_eq!(graph.reverse_cuthill_mckee(), [4, 3, 2, 1, 0]);

// The induced subgraph on {2, 1}, renumbered in the given order.
let subgraph = graph.subgraph([2, 1]);
assert_eq!(subgraph.row(0), [1]); // old vertex 2 is new vertex 0
assert_eq!(subgraph.row(1), [0]); // old vertex 1 is new vertex 1
```

`connected_components` requires a square graph (`n() == m()`) and numbers the
components in discovery order: the component of vertex `0` is `0`, then the first
unvisited vertex starts component `1`, and so on.

## Reordering

Reordering changes vertex, element, and face numbering without changing the geometry or
the connectivity. It is used to improve memory locality and to reduce the bandwidth of
the element graph.

The bandwidth of an element connectivity is the largest difference between the maximum
and minimum vertex index of an element:

$$
b = \max_{K} \left( \max_{v \in K} v - \min_{v \in K} v \right).
$$

`bandwidth(elems)` returns the maximum bandwidth and the mean bandwidth over the
elements.

| Method | Ordering |
| --- | --- |
| `reorder_rcm()` | reverse Cuthill-McKee on the vertex-to-vertex graph |
| `reorder_hilbert()` | Hilbert curve over the vertex coordinates |
| `random_shuffle()` | uniformly shuffled vertices, elements, and faces (fixed seed) |

The Hilbert ordering normalizes the vertex coordinates by the bounding box, quantizes
each coordinate on a $2^{16}$ grid, and linearizes the result with the Hilbert curve of
the `lindel` library. It is a geometric ordering: spatially close vertices receive
close indices, which gives good cache locality, but it does not minimize the bandwidth
as well as RCM.

`reorder_rcm` and `reorder_hilbert` return a tuple
`(mesh, vert_ids, elem_ids, face_ids)`. The three arrays map **new to old** indices:
`new_mesh.vert(k) == old_mesh.vert(vert_ids[k])`, and similarly for elements and faces.
The new element and face orders are obtained by sorting on the minimum vertex index.

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, bandwidth, rectangle_mesh};

let msh = rectangle_mesh::<Mesh2d>(1.0, 10, 2.0, 15);

// Shuffling makes the initial bandwidth large.
let shuffled = msh.random_shuffle();
let (b_max, _b_mean) = bandwidth(shuffled.elems());

// RCM minimizes the bandwidth of the element connectivity.
let (rcm, vert_ids, elem_ids, face_ids) = shuffled.reorder_rcm();
let (b_rcm, _b_mean_rcm) = bandwidth(rcm.elems());

assert!(b_rcm < b_max);
assert_eq!(vert_ids.len(), msh.n_verts());
assert_eq!(elem_ids.len(), msh.n_elems());
assert_eq!(face_ids.len(), msh.n_faces());
```

`random_shuffle` uses the fixed seed `1234`, so a given mesh always produces the same
shuffled mesh. This is what makes reordering benchmarks reproducible.

## Partitioning

The element partition is computed from the element-to-element graph, which is built
from the internal faces (`element_pairs(all_faces)`). It is exposed by the
`Partitioner` trait and by the `Mesh::partition` method:

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, partition::HilbertPartitioner, rectangle_mesh};

let mut msh = rectangle_mesh::<Mesh2d>(1.0, 10, 2.0, 15);

let (quality, imbalance) = msh
	.partition::<HilbertPartitioner>(4, None)
	.unwrap();

assert!(quality < 0.2);
assert!(imbalance <= 1e-12);
```

`partition(n_parts, weights)` returns two diagnostics and writes the partition into the
element tags: elements of part `i` receive the tag `i + 1`.

| Partitioner | Cargo feature | Algorithm |
| --- | --- | --- |
| `HilbertPartitioner` | none | sort elements by the Hilbert index of their centers, then cut into chunks |
| `RCMPartitioner` | none | RCM ordering of the element graph, then cut into chunks |
| `KMeansPartitioner2d` | `coupe` | Hilbert initialization followed by weighted k-means on element centers |
| `KMeansPartitioner3d` | `coupe` | same as above in 3D |
| `MetisPartitioner<MetisRecursive>` | `metis` | recursive bisection with METIS |
| `MetisPartitioner<MetisKWay>` | `metis` | k-way partitioning with METIS |

For the geometric and RCM partitioners, the elements are first ordered (by Hilbert
index or by RCM), then cut into contiguous chunks of similar weight. If the cumulative
weight before element $j$ is $c_j$, the total weight is $W$, and the weight of element
$j$ is $w_j$, the part index is

$$
p_j = \min\left( \left\lfloor n_{\text{parts}}\,
\frac{c_j + w_j/2}{W} \right\rfloor,\; n_{\text{parts}} - 1 \right).
$$

With uniform weights this produces $n_{\text{parts}}$ chunks whose sizes differ by at
most one element. The `coupe` and `metis` partitioners instead work directly on the
element graph or on the element centers.

### Quality and imbalance

Both diagnostics are returned by `partition`:

- the **quality** is the fraction of element-adjacency entries that connect different
  parts:

$$
q = \frac{\#\{(i,j) : i \sim j,\ p_i \neq p_j\}}{\#\{(i,j) : i \sim j\}}.
$$

  It is $0$ when no internal face is cut and grows towards $1$ when the parts are
  interleaved. Since the adjacency graph is symmetric, every cut face is counted twice,
  which does not change the ratio.
- the **imbalance** is the normalized spread of the part weights:

$$
\mathrm{imb} = \frac{\max_i W_i - \min_i W_i}{\frac{1}{n_{\text{parts}}} \sum_i W_i},
$$

  where $W_i$ is the sum of the element weights in part $i$. It is $0$ for perfectly
  balanced parts.

With uniform weights on a 10 by 15 rectangle split into four parts, Hilbert gives
$q \approx 0.0986$ and $\mathrm{imb} = 0$; RCM gives $q \approx 0.0930$ and
$\mathrm{imb} = 0$ on the same mesh.

### Extracting a part

`get_partition(i)` extracts part `i` (tag `i + 1`) as a `SubMesh`:

```rust
# extern crate tmesh;
use tmesh::mesh::{Mesh, Mesh2d, partition::HilbertPartitioner, rectangle_mesh};

let mut msh = rectangle_mesh::<Mesh2d>(1.0, 10, 2.0, 15);
msh.partition::<HilbertPartitioner>(4, None).unwrap();

let part0 = msh.get_partition(0);
assert_eq!(part0.mesh.n_elems(), 63);
// The arrays map the submesh entities back to the parent mesh.
assert_eq!(part0.parent_elem_ids.len(), part0.mesh.n_elems());
assert_eq!(part0.parent_vert_ids.len(), part0.mesh.n_verts());
assert_eq!(part0.parent_face_ids.len(), part0.mesh.n_faces());
```

`SubMesh::new(&msh, filter)` is the general extraction routine: it keeps the elements
whose tag satisfies `filter`, together with their vertices and faces, and records the
parent indices.

> [!WARNING]
> `partition` and `set_partition` **overwrite the element tags** with the part index
> plus one. Any region tag used to mark materials or subdomains is lost. Save the tags
> first if they are needed afterwards.

`pytucanos` exposes the same partitioners:

```python
import numpy as np
from pytucanos import Mesh2d, PartitionerType

msh = Mesh2d.rectangle_mesh(np.linspace(0.0, 1.0, 10), np.linspace(0.0, 2.0, 15))
quality, imbalance = msh.partition(4, PartitionerType.Hilbert)

# Elements of part i are tagged i + 1.
assert np.unique(msh.get_etags()).tolist() == [1, 2, 3, 4]
```

> [!TIP]
> **Graph and partitioner error and panic messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `The number of partitions must be > 0` | `n_parts` is zero. |
> | `Invalid number of weights: {len} (expected {n_elems})` | The optional weight array must provide one weight per element. |
> | `The weights must be finite and >= 0` | A weight is `NaN`, infinite, or negative. |
> | `Partitioner only available for D=2` / `Partitioner only available for D=3` | `KMeansPartitioner2d` and `KMeansPartitioner3d` are dimension-specific; use the matching one. |
> | `too many connected components` | `connected_components` overflowed its component counter. This cannot happen for a realistic mesh. |
> | METIS error string | The METIS backend rejected the graph. Check the element graph and the requested number of parts. |
>
> Panics to be aware of:
>
> | Panic | Cause |
> | --- | --- |
> | ``called `Option::unwrap()` on a `None` value`` | `CSRGraph::new` was called with an empty `indices` array. |
> | ``assertion `left != right` failed`` | `from_edges` received a duplicate edge or a self-loop, which would create a duplicate column. |
> | ``assertion `left == right` failed`` in `connected_components` | The graph is not square (`n() != m()`). |
> | ``assertion `left == right` failed`` in `set_partition` | The partition array does not have one entry per element. |
> | ``called `Option::unwrap()` on a `None` value`` in `row_and_values` | The graph has no attached values; use `row` instead. |

`q` is `NaN` when the mesh has no internal face (a single element), because the
denominator of the quality formula is zero.
