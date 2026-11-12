# Geometry

This page describes how `tucanos` links mesh entities to a geometric model and how geometric quantities are computed.

## Notation

- spatial dimension: $d$,
- mesh vertex: $\mathbf{x}_i \in \mathbb{R}^d$,
- topological tag: $\tau=(q,t)$ with topological dimension $q$ and integer tag $t$.

Projection notation used on this page:

- tag-routed projection map: $\Pi_{\tau}(\mathbf{x})$,
- projection distance: $\operatorname{dist}(\mathbf{x},\tau)=\|\mathbf{x}-\Pi_{\tau}(\mathbf{x})\|$.

The geometric model is queried with a point and a topological tag to select the relevant geometric entity (surface, curve, or corner).

## Problem solved by geometry handling

Given a mesh and topological labels, geometry handling provides:

- projection of mesh points to the geometric support,
- geometric error indicators (distance and normal-angle mismatch),
- curvature directions on tagged patches,
- quadratic elevation by inserting new points and projecting them.

## Inputs and outputs

Inputs:

- a mesh with consistent tags,
- a topology map (`MeshTopology`) that assigns topological entities/tags,
- a geometric representation (none, analytical, or meshed boundary-based).

Outputs:

- projected vertices,
- projection distances,
- normal-angle deviations,
- optional curvature vectors per patch.

## Projection semantics

For a point $\mathbf{x}$ and topological tag $\tau$:

1. Select the geometric entity associated with $\tau$.
2. Compute the closest projection $\Pi_{\tau}(\mathbf{x})$.
3. Return the projection distance $\operatorname{dist}(\mathbf{x},\tau)$.

The same dispatch rule is used during bulk projection of boundary vertices.

## Geometric quality indicators

Two core indicators are used:

- distance error: maximum projection distance over selected mesh points,
- normal-angle error: maximum angle between discrete mesh normal and geometry normal.

These indicators are used to validate geometric fidelity after projection and after quadratic elevation.

## Meshed geometry representation

For meshed geometries, `tucanos` stores boundary patches by topological tags and uses spatial indices for nearest-point queries.

High-level algorithm:

1. Build tagged patch subsets from a boundary mesh.
2. Build per-patch spatial indices.
3. At query time, route by topological tag and project on the corresponding patch index.

This enables robust projection on piecewise linear geometric supports.

## Quadratic elevation

Quadratic conversion is geometry-aware:

1. Insert mid-edge nodes from linear connectivity.
2. Classify inserted nodes by topological support.
3. Project boundary-support nodes on geometry.
4. Keep interior nodes unchanged.

Output: a quadratic mesh whose curved boundary better matches the geometric model.

## Rust example

```rust
use tmesh::mesh::{Mesh, Mesh2d, rectangle_mesh};
use tucanos::mesh::topology::MeshTopology;

let mut msh = rectangle_mesh::<Mesh2d>(1.0, 10, 2.0, 15);
msh.fix().unwrap();
let topo = MeshTopology::new(&msh);

// Geometry objects use topo tags to project/query boundary entities.
// Typical outputs are projection distance and max distance/angle indicators.
let _ = topo;
```

## Python note

In `pytucanos`, geometry-aware operations follow the same logic: create/fix the mesh, build or provide the geometry model, then evaluate distances/angles or run projection-dependent workflows.