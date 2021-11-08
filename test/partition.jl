using Test

using GCMMesh
using GCMMesh.Mesh
using GCMMesh.Partition


FT = Float64

@testset "test rsb_partition" begin
    x1min, x1max = FT(0), FT(1)
    x2min, x2max = FT(0), FT(1)
    limits = (x1min, x2min, x1max, x2max)
    n1, n2 = 4, 4 # # of elements along x1 and x2 axis respectively
    npart = 5 # # of partitions

    mesh = equispaced_rectangular_mesh(limits..., n1, n2, (false, false))
    partition = rsb_partition(mesh, npart)

    @test extrema(partition) == (0, npart - 1) # verify if process numbers are assigned to all elems
    partlen = zeros(Int, npart)
    for i = 1:npart
        partlen[i] = count(partition .== i - 1)
    end
    @test sum(partlen) == mesh.nelems             # verify sum(nlocalelems) == nelems 
    @test maximum(partlen) - minimum(partlen) ≤ 1 # check if partion lengths are even
end
