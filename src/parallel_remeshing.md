# Parallel remeshing

At each level, the parallel remesher:

1. Partitions the elements and builds one submesh per partition. Each partition starts from a clone of the global topology. Only partition-interface face tags are cleared with `Topology::clear(filter)` before `MeshTopology::new_from` rebuilds the local hierarchy; see [Reusing an existing hierarchy](topology.md#reusing-an-existing-hierarchy-new_from).
2. Remeshes each partition independently with the same local parameters.
3. Extracts the partition interfaces and remeshes them recursively while the interface is large enough, or with the sequential remesher at the finest level.
4. Merges the remeshed partitions and interfaces, then removes the temporary interface boundary tags.

Remeshing failures are wrapped as `Failed to remesh partition {i_part}.\n{error}`; see [Validity](validity.md#remesher-consistency-errors). Partitioner input errors are `The number of partitions must be > 0`, `Invalid number of weights: {n} (expected {n_elems})`, and `The weights must be finite and >= 0`.
