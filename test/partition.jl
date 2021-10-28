using Test

using GCMMesh
using GCMMesh.Mesh
using GCMMesh.Partition


FT = Float64

x1min, x1max = FT(0), FT(1)
x2min, x2max = FT(0), FT(1)
limits = (x1min, x2min, x1max, x2max)
n1, n2 = 4, 4 # # of elements along x1 and x2 axis respectively
npart = 5 # # of partitions

mesh = equispaced_rectangular_mesh(limits..., n1, n2, (false, false))

partition = rsb_partition(mesh, npart)
