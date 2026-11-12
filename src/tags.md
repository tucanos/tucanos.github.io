# Tags

## Element and face tags

Mesh entities are tagged as:

- element tags $T_{K,i}$ on simplices $K_i$,
- face tags $T_{F,i}$ on boundary/interface faces $F_i$.

Face tags are part of mesh validity, not optional annotations.

Tag values follow these conventions:

- tags are signed integers,
- value $0$ is avoided (automatic tag creation uses nonzero values),
- in remeshing workflows, a negative face tag means this face is protected and should not be modified.

## Face tagging rules

A face must be tagged when:

- it belongs to a single element (external boundary),
- it separates two elements with different element tags,
- it is shared by more than two elements (non-manifold interface).

A face shared by exactly two elements with the same element tag must not be tagged.

## Face tagging errors

Face-tagging rules are enforced by `Mesh::check()`.

> [!TIP]
> **`Mesh::check()` messages**
>
> | Message | Meaning / action |
> | --- | --- |
> | `Boundary face {f:?} not tagged` | A face adjacent to a single element is missing from the face list. Call `fix()` or add it. |
> | `Internal boundary face {f:?} not tagged (elems = ..., tags = ...)` | A face shared by elements with different tags (or by more than two elements) is untagged. Tag it, or make element tags consistent. |
> | `Tagged face inside the domain: center = {fc:?}` | A face shared by exactly two elements with the same element tag is tagged. Remove it from the face list, or change the element tags so the interface is real. `fix()` does not remove it. |

## Transferring tags

`transfer_tags(src, dst)` transfers element tags when both meshes have the same dimension, or face tags when the destination has dimension `src_dim - 1`.

> [!TIP]
> **Tag transfer error**
>
> Any other dimension combination panics with `Invalid mesh element dimension`. Transfer element tags between meshes of equal dimension, or transfer element tags to face tags from a mesh of dimension $d$ to one of dimension $d-1$.
