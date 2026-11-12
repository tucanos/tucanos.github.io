# Non-manifold meshes

`tmesh` accepts non-manifold simplex meshes when non-manifold interfaces are explicitly represented and tagged.

## Detection rule

In face-to-element connectivity, a face is non-manifold when it has more than two adjacent elements.

## Required handling

Non-manifold faces are treated as mandatory internal boundaries:

- they must appear in the face list,
- they must have face tags,
- missing non-manifold interface faces/tags may be inserted by `fix()`.

Mesh checks fail if such faces are untagged or inconsistent.

## Tagging consistency

Face-tag policy for all interface types is:

- 1 adjacent element: boundary face, tag required,
- 2 adjacent elements with different element tags: interface face, tag required,
- 2 adjacent elements with same element tag: no interface tag,
- more than 2 adjacent elements: non-manifold interface, tag required.

## Current limitation

For some `ids.len() > 2` auto-tagging paths with repeated adjacent element tags, the implementation currently contains an unimplemented branch.

Practical consequence: for complex non-manifold configurations, prepare face connectivity/tagging explicitly before final validation.
