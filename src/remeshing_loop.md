# The remeshing algorithm

The remesher applies local cavity-based operators to move the current mesh toward a quasi-unit mesh with respect to the prescribed metric field.

Implementation reference: `tucanos/src/remesher/sequential.rs` (`Remesher::remesh`).

## Scheduler and operation order

The sequential scheduler executes an ordered list of steps:

1. collapse,
2. split,
3. swap,
4. smooth.

In the default parameter builder (`RemesherParams::new`), this pattern is repeated with two swap passes (different quality targets), then finalized by two additional swap passes. In `new_ordered`, the pattern is collapse → split → ordered swap → smooth.

This ordering is important: collapse removes clearly over-refined edges first, split restores under-resolved areas, swap improves connectivity/quality, and smoothing relaxes the vertex placement while preserving validity.

## Why topology hierarchy is central

All operators use topology tags and parent relations (from `MeshTopology`) to prevent invalid geometric moves:

- edge tags are inferred from vertex tags, tagged faces, and parent sets (`edge_tag`),
- candidate reconstructions are rejected when parent/child compatibility is violated (`topo.parent(...)` checks),
- additional cavity-level topological checks prevent element/face pinching and duplicate entities (`FilledCavity::check_topo`).

Without this hierarchy, operations that look geometrically valid can still break feature consistency at ridges/corners/interfaces.

## Collapse step

Implementation reference: `tucanos/src/remesher/collapse.rs` (`collapse`, `perform_collapse`).

### Loop order

Edges are sorted by:

1. topological dimension of edge tag,
2. increasing metric edge length.

Then each edge is tested in this order.

### Cavity creation

For candidate edge $(i_0,i_1)$, collapse is evaluated from a **vertex cavity** built around `i0` (`Cavity::init_from_vertex`).

The cavity is then filled with existing vertex `i1` (`FilledCavityType::ExistingVertex`).

### Acceptance conditions

An edge collapse is accepted only if all checks pass:

1. edge length below threshold `params.l`,
2. edge/tag compatibility with topology hierarchy,
3. no fixed-boundary constraint violation (negative tags),
4. cavity reconstruction passes topological checks,
5. boundary normals remain geometry-consistent (`check_normals`),
6. new elements satisfy quality/length bounds (`FilledCavity::check`).

If accepted, all old cavity elements are deleted, vertex `i0` is removed, and reconstructed elements/tagged faces are inserted.

## Split step

Implementation reference: `tucanos/src/remesher/split.rs` (`split`, `try_split`, `extend_cavity`, `edge_split`).

### Loop order

Edges are sorted by:

1. parent-topology dimension,
2. decreasing metric edge length.

Only edges above the split threshold are processed.

### Cavity creation

For candidate edge $(i,j)$:

1. build an **edge cavity** (`Cavity::init_from_edge`),
2. create candidate split point from seed barycenter,
3. project this point to geometry when edge lies on boundary/interface,
4. build filled cavity with `FilledCavityType::EdgeCenter`.

### Acceptance conditions

A split is accepted only if:

1. edge is not fixed (`tag < 0` is rejected),
2. cavity passes quality/length bounds,
3. new boundary normals stay within allowed angle.

### Cavity expansion and local deletion fallback

When the direct split is invalid on boundary configurations:

- a cavity expansion loop (`extend_cavity`) attempts to grow the cavity through invalid faces up to `max_extensions`,
- in specific 3D boundary cases, one boundary element can be removed (`remove_boundary_element`) before retrying.

If accepted, original cavity elements are deleted and replaced by split reconstruction.

## Edge swap step

Implementation reference: `tucanos/src/remesher/swap.rs` (`swap`, `find_swap_vertex`, `perform_swap`).

### Loop order

All edges are iterated (unordered mode), or with an ordered strategy when `SwapParams::ordered` is enabled.

### Cavity creation

For each edge:

1. build an **edge cavity** (`Cavity::init_from_edge`),
2. enumerate cavity vertices (except edge endpoints) as possible replacement vertices,
3. build candidate filled cavity with `FilledCavityType::ExistingVertex(candidate)`.

### Acceptance conditions

A swap candidate is accepted only if:

1. local cavity quality is below target (otherwise no-op),
2. edge is not fixed and not forbidden by tag dimension,
3. topology parent checks are valid for candidate vertex,
4. cavity topological checks pass,
5. geometry-normal consistency check passes,
6. quality/length checks pass,
7. among valid candidates, it maximizes the new minimum cavity quality.

If accepted, old cavity elements are deleted and replaced by swapped reconstruction.

> [!WARNING]
> Face-swap is not implemented in the current remeshing loop. Only edge-based swap is available.

## Smoothing step

Implementation reference: `tucanos/src/remesher/smooth.rs` (`smooth`, `smooth_iter`).

### Loop order

Smoothing runs for `n_iter` passes over vertices.

### Cavity creation

For each vertex:

1. build a **vertex cavity** (`Cavity::init_from_vertex`),
2. compute admissible neighbor set using topology hierarchy (`get_smoothing_neighbors`),
3. compute candidate smoothed position (Laplacian / Laplacian2 / Avro),
4. apply relaxation factors and project to geometry when needed,
5. build candidate filled cavity with `FilledCavityType::MovedVertex`.

### Acceptance conditions

A smoothing move is accepted only if:

1. vertex is not fixed (`dim=0` or negative tag are skipped),
2. candidate respects geometry-normal angle bound,
3. candidate does not reduce minimum cavity quality.

If accepted, vertex position is updated and metric at moved vertex is re-interpolated from cavity elements (barycentric interpolation in best containing/closest element).
