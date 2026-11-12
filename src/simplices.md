# Simplices


## Topological versus geometric simplex

`tmesh` uses two complementary abstractions:

- `Simplex`: topological simplex (indices, local numbering, incidence),
- `GSimplex<D>`: geometric simplex in an ambient space of dimension $D$ (coordinates and geometric operators).

The link between both is explicit in the associated types:

- for a topological simplex `C: Simplex`, `C::GEOM<D>` is its geometric representation,
- for a geometric simplex `G: GSimplex<D>`, `G::TOPO` is its topological representation.

`GSimplex` provides geometric operations on element coordinates, including:

- integration (`integrate`),
- measure (`vol`),
- normals when defined for the embedding (`normal`),
- characteristic radius (`radius`),
- barycentric mapping (`bcoords`, `vert`),
- projection and distance (`project`, `distance`),
- quality and geometric bounds (`gamma`, `bounding_box`).

This separation is what allows, for example, a triangle topology to be embedded in different ambient spaces while keeping the same topological conventions.

## Core assumptions from `Simplex`

A simplex is represented by an ordered tuple of vertex indices:

- edge: $(v_0,v_1)$,
- triangle: $(v_0,v_1,v_2)$,
- tetrahedron: $(v_0,v_1,v_2,v_3)$.

The order is not arbitrary: it defines orientation and local numbering.

The key face convention is:

- face #$i$ is the face opposite to vertex #$i$,
- equivalently, face #$i$ does not contain vertex #$i$.

This is enforced both topologically (`Simplex::face(i)`) and geometrically (`GSimplex::face(i)`).

## Index storage assumptions (`Idx`)

Simplex indices are stored through the `Idx` abstraction. In `tmesh`, supported index types are:

- unsigned 32-bit: `u32`,
- signed 32-bit: `i32`,
- signed 64-bit: `i64`,
- platform-native unsigned: `usize`.

So indices can be stored using signed or unsigned integers, with 32-bit support and 64-bit support through signed `i64` (and `usize` on 64-bit platforms).

In all cases, indices must be convertible to/from `usize` for internal connectivity operations.

## Edge local numbering

For an edge $(v_0,v_1)$:

- the unique edge entity is $(v_0,v_1)$,
- faces are 0D nodes (endpoints),
- opposite-vertex convention gives:
	- face 0: $(v_1)$,
	- face 1: $(v_0)$.

So face $i$ excludes vertex $v_i$ exactly, just as for higher-order simplices.

Orientation is given by endpoint order; swapping endpoints flips orientation.

## Orientation assumptions

Orientation is encoded by vertex order.

- For edges, orientation follows endpoint order.
- For triangles and tetrahedra, local vertex ordering defines the orientation used by incidence and geometric operators.

Swapping two vertices flips orientation. In the trait implementation, `invert()` swaps local vertices 0 and 1 to reverse orientation.

In mesh workflows, valid elements are expected to have positive measure after orientation fixes.

## Triangle local numbering

For a triangle $(v_0,v_1,v_2)$:

- edges: $(v_0,v_1)$, $(v_1,v_2)$, $(v_2,v_0)$,
- faces (edge-faces) by opposite-vertex convention:
	- face 0: $(v_1,v_2)$,
	- face 1: $(v_2,v_0)$,
	- face 2: $(v_0,v_1)$.

So face $i$ excludes vertex $v_i$ exactly.

## Tetrahedron local numbering

For a tetrahedron $(v_0,v_1,v_2,v_3)$:

- edges:
	- $(v_0,v_1)$, $(v_1,v_2)$, $(v_2,v_0)$,
	- $(v_0,v_3)$, $(v_1,v_3)$, $(v_2,v_3)$,
- triangular faces by opposite-vertex convention:
	- face 0: $(v_1,v_2,v_3)$,
	- face 1: $(v_2,v_0,v_3)$,
	- face 2: $(v_0,v_1,v_3)$,
	- face 3: $(v_0,v_2,v_1)$.

The vertex order inside each face is chosen consistently with the parent-tetra orientation convention used in `tmesh`.

## Equality and canonicalization assumptions

Two different notions are used:

- `is_same`: orientation-aware simplex equality under valid local permutations,
- `sorted`: canonical index ordering for hashing/lookup (orientation is not preserved by sorting).

In practice:

- use orientation-preserving checks when geometric sign matters,
- use sorted forms for connectivity maps and deduplication keys.

