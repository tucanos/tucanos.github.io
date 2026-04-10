# Topology

This module provides data structures and algorithms to represent, build, and query the topological relationships of a mesh. It extracts hierarchical adjacency information from mesh elements (e.g., tetrahedrons, triangles) and faces, mapping out how entities of different dimensions (volumes, faces, edges, vertices) connect to one another.

## Core Concepts

* **`Dim`**: Represents the dimension of an entity (e.g., `0` for vertices, `1` for edges, `2` for faces, `3` for volumes).
* **`Tag`**: A unique integer identifier for a specific entity within a given dimension.
* **`TopoTag`**: A tuple of `(Dim, Tag)` uniquely identifying any topological entity across the entire mesh.
* **Hierarchy**: In this system, "parents" are strictly higher-dimensional entities that bound a "child" entity. For example, a 3D Volume (Dim 3) is a parent to its bounding 2D Faces (Dim 2). A 2D Face is a parent to its 1D Edges.

---

## Core Data Structures

### `TopoNode`
Represents a single topological entity in the mesh and its immediate relationships.
```rust
pub struct TopoNode {
    pub tag: TopoTag,            // The unique identifier (Dim, Tag)
    pub children: HashSet<Tag>,  // Tags of lower-dimensional boundary entities (Dim - 1)
    pub parents: HashSet<Tag>,   // Tags of higher-dimensional bounding entities (Dim + 1)
}
```

### `Topology`
The central registry containing all `TopoNode`s, categorized by dimension, and caching their hierarchical relationships.
```rust
pub struct Topology {
    dim: Dim,
    entities: Vec<Vec<TopoNode>>, 
    parents: FxHashMap<(TopoTag, TopoTag), TopoTag>,
}
```
* **`entities`**: A 2D vector where `entities[d]` holds all the `TopoNode`s of dimension `d`.
* **`parents`**: A cached mapping used to quickly find the closest common parent of a given pair of `TopoTag`s.

### `MeshTopology`
A high-level wrapper that binds a `Topology` to a specific `Mesh` instance. It manages the `Topology` and stores the resolved topological tags for all vertices (`vtags`).
```rust
pub struct MeshTopology {
    topo: Topology,
    vtags: Vec<TopoTag>,
}
```

---

## API Reference: `Topology`

### Initialization and Serialization
* **`new(dim: Dim) -> Self`**: Initializes an empty topology supporting up to the specified dimension.
* **`from_json(fname: &str) -> Result<Self>`**: Deserializes a topology from a JSON file and automatically recomputes the parent caches.
* **`to_json(&self, fname: &str) -> Result<()>`**: Serializes the topology to a pretty-printed JSON file.

### Querying Tags and Entities
* **`ntags(&self, dim: Dim) -> usize`**: Returns the number of tags registered in a specific dimension.
* **`tags(&self, dim: Dim) -> Vec<Tag>`**: Returns a list of all tags present in a specific dimension.
* **`first_available_tag(&self, dim: Dim) -> Tag`**: Calculates the next available unique tag for a given dimension (the maximum absolute tag value + 1).
* **`get(&self, tag: TopoTag) -> Option<&TopoNode>`**: Retrieves a reference to a `TopoNode` by its `TopoTag`.

### Hierarchy Construction
* **`insert(&mut self, tag: TopoTag, parents: &[Tag])`**: Manually inserts a new entity. It automatically updates the `children` set of the specified parent entities (dimension `tag.0 + 1`).
* **`get_from_parents(&self, dim: Dim, parents: &[Tag]) -> Option<&TopoNode>`**: Finds a node in the specified dimension that has the *exact* set of given parent tags.
* **`parent(&self, topo0: TopoTag, topo1: TopoTag) -> Option<TopoTag>`**: Retrieves the closest common parent of two topology entities, leveraging the internal `FxHashMap` cache.

### Automatic Topology Inference
The module contains robust algorithms to automatically deduce lower-dimensional topology from a mesh's elements and boundary faces:
* **`update_from_elems_and_faces`**: The primary engine used by `MeshTopology`. It iterates over the mesh elements and faces, infers the tags for intermediate dimensions (like edges in a 3D mesh), resolves tagging conflicts at corners/junctions (`check_and_fix`), and returns a fully resolved vector of vertex tags (`vtags`).
* **`clear<F>(&mut self, filter: F)`**: Filters out topology nodes based on a closure, cleaning up dangling child/parent relationships along the way.

---

## API Reference: `MeshTopology`

* **`new<const D: usize, M: Mesh<D>>(msh: &M) -> Self`**: Constructs a full topological mapping from a generic mesh `M`. It analyzes element tags (`etags`) and face tags (`ftags`), generating unique tags for unmarked internal boundaries, edges, and vertices.
* **`topo(&self) -> &Topology`**: Grants immutable access to the underlying `Topology`.
* **`vtags(&self) -> &[TopoTag]`**: Returns a slice containing the computed `TopoTag` for every vertex in the underlying mesh, indexed by the vertex ID.

---

## Usage Example

While usually constructed automatically via `MeshTopology::new(&mesh)`, the topology can be constructed manually (which is often useful for testing and custom mesh generation):

```rust
let mut t = Topology::new(3);

// Insert a 3D volume (Dim 3, Tag 1) with no parents
t.insert((3, 1), &[]);

// Insert 2D faces (Dim 2) whose parent is the volume (Tag 1)
t.insert((2, 1), &[1]);
t.insert((2, 2), &[1]);

// Insert a 1D edge (Dim 1) bounded by both faces
t.insert((1, 1), &[1, 2]); // Edge shared by Face 1 and Face 2

// Automatically resolve all parent-child hierarchy caches
t.compute_parents();

// Fetch an edge bounded by faces 1 and 2
let edge = t.get_from_parents(1, &[1, 2]).unwrap();
assert_eq!(edge.tag.1, 1);
```