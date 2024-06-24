# Remesher

Input metric                 |  Adapted mesh
:---------------------------:|:-------------------------:
![metric](images/metric.svg) |  ![adapted mesh](images/adapted.svg)

The idea is to perform a series of local modification on an input mesh to have a 
- edge lengths that are of approximately unit length in the metric space
- element qualities as high as possible

In order to ensure the validity of the final mesh, validity is enforced at each step.
