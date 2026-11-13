# Parallel remeshing

The parallel remesher partitions the mesh, remeshes each partition independently, and
stitches the partitions back together. Partition-selection and remeshing failures are
reported with the messages listed in [Validity](validity.md#remesher-consistency-errors).
