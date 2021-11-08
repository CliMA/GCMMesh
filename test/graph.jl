using Test
using Plots
using LinearAlgebra

using GCMMesh
using GCMMesh.Mesh

using GCMMesh.Graphs

FT = Float64

x1min, x1max = FT(0), FT(1)
x2min, x2max = FT(0), FT(1)

limits = (x1min, x2min, x1max, x2max)

function is_in_fiedler_vector_space(
    eigvec,
    graph::Graph{I,FT},
    evec_tol = FT(1.0E-6),
) where {I,FT}
    nv = graph.nverts
    vl, vu = FT(1), FT(1) # not used
    il, iu = 1, nv
    st, en = 2, 2
    abstol = FT(1.0E-14)    # abs tolerance for syevr! LAPACK solver
    eig_reltol = FT(1.0E-4) # relative tolerance for finding fiedler eigenvalues with algebraic
                            # multiplicity > 1
    lmat = Graphs.laplacian_matrix_full(graph)
    eigvals, eigvecs = LinearAlgebra.LAPACK.syevr!('V', 'I', 'L', lmat, vl, vu, il, iu, abstol)
    
    fdlr_eig = eigvals[2]
    for i in 3:nv
        if abs(fdlr_eig - eigvals[i]) < eig_reltol
            en = i
        else
            break
        end
    end
    diff = zeros(FT, nv)
    for i in st:en
        diff .+= (eigvec' * eigvecs[:, i]) .* eigvecs[:, i]
    end
    return norm(diff .- eigvec) < evec_tol
end

@testset "graph for a simple 4 x 4 rectangular grid" begin
    @testset "face connectivity graph for a 4×4 element quad mesh" begin
        mesh = equispaced_rectangular_mesh(limits..., 4, 4, (false, false))
        fgraph = build_face_graph(mesh)
        off = fgraph.edge_offset

        @test fgraph.nverts == 4 * 4
        @test fgraph.edge_data[off[1]:off[2]-1] == [2, 5]
        @test fgraph.edge_data[off[2]:off[3]-1] == [1, 3, 6]
        @test fgraph.edge_data[off[3]:off[4]-1] == [2, 4, 7]
        @test fgraph.edge_data[off[4]:off[5]-1] == [3, 8]
        @test fgraph.edge_data[off[5]:off[6]-1] == [1, 6, 9]
        @test fgraph.edge_data[off[6]:off[7]-1] == [2, 5, 7, 10]
        @test fgraph.edge_data[off[7]:off[8]-1] == [3, 6, 8, 11]
        @test fgraph.edge_data[off[8]:off[9]-1] == [4, 7, 12]
        @test fgraph.edge_data[off[9]:off[10]-1] == [5, 10, 13]
        @test fgraph.edge_data[off[10]:off[11]-1] == [6, 9, 11, 14]
        @test fgraph.edge_data[off[11]:off[12]-1] == [7, 10, 12, 15]
        @test fgraph.edge_data[off[12]:off[13]-1] == [8, 11, 16]
        @test fgraph.edge_data[off[13]:off[14]-1] == [9, 14]
        @test fgraph.edge_data[off[14]:off[15]-1] == [10, 13, 15]
        @test fgraph.edge_data[off[15]:off[16]-1] == [11, 14, 16]
        @test fgraph.edge_data[off[16]:off[17]-1] == [12, 15]
    end

    @testset "face connectivity graph for a 4×4 element quad mesh, periodic boundaries" begin
        mesh = equispaced_rectangular_mesh(limits..., 4, 4, (true, true))
        fgraph = build_face_graph(mesh)
        off = fgraph.edge_offset

        @test fgraph.nverts == 4 * 4
        @test fgraph.edge_data[off[1]:off[2]-1] == [2, 4, 5, 13]
        @test fgraph.edge_data[off[2]:off[3]-1] == [1, 3, 6, 14]
        @test fgraph.edge_data[off[3]:off[4]-1] == [2, 4, 7, 15]
        @test fgraph.edge_data[off[4]:off[5]-1] == [1, 3, 8, 16]
        @test fgraph.edge_data[off[5]:off[6]-1] == [1, 6, 8, 9]
        @test fgraph.edge_data[off[6]:off[7]-1] == [2, 5, 7, 10]
        @test fgraph.edge_data[off[7]:off[8]-1] == [3, 6, 8, 11]
        @test fgraph.edge_data[off[8]:off[9]-1] == [4, 5, 7, 12]
        @test fgraph.edge_data[off[9]:off[10]-1] == [5, 10, 12, 13]
        @test fgraph.edge_data[off[10]:off[11]-1] == [6, 9, 11, 14]
        @test fgraph.edge_data[off[11]:off[12]-1] == [7, 10, 12, 15]
        @test fgraph.edge_data[off[12]:off[13]-1] == [8, 9, 11, 16]
        @test fgraph.edge_data[off[13]:off[14]-1] == [1, 9, 14, 16]
        @test fgraph.edge_data[off[14]:off[15]-1] == [2, 10, 13, 15]
        @test fgraph.edge_data[off[15]:off[16]-1] == [3, 11, 14, 16]
        @test fgraph.edge_data[off[16]:off[17]-1] == [4, 12, 13, 15]
    end

    @testset "vertex connectivity graph for a 4×4 element quad mesh" begin
        mesh = equispaced_rectangular_mesh(limits..., 4, 4, (false, false))
        vgraph = build_vertex_graph(mesh)
        off = vgraph.edge_offset

        @test vgraph.nverts == 4 * 4
        @test vgraph.edge_data[off[1]:off[2]-1] == [2, 5, 6]
        @test vgraph.edge_data[off[2]:off[3]-1] == [1, 3, 5, 6, 7]
        @test vgraph.edge_data[off[3]:off[4]-1] == [2, 4, 6, 7, 8]
        @test vgraph.edge_data[off[4]:off[5]-1] == [3, 7, 8]
        @test vgraph.edge_data[off[5]:off[6]-1] == [1, 2, 6, 9, 10]
        @test vgraph.edge_data[off[6]:off[7]-1] == [1, 2, 3, 5, 7, 9, 10, 11]
        @test vgraph.edge_data[off[7]:off[8]-1] == [2, 3, 4, 6, 8, 10, 11, 12]
        @test vgraph.edge_data[off[8]:off[9]-1] == [3, 4, 7, 11, 12]
        @test vgraph.edge_data[off[9]:off[10]-1] == [5, 6, 10, 13, 14]
        @test vgraph.edge_data[off[10]:off[11]-1] == [5, 6, 7, 9, 11, 13, 14, 15]
        @test vgraph.edge_data[off[11]:off[12]-1] == [6, 7, 8, 10, 12, 14, 15, 16]
        @test vgraph.edge_data[off[12]:off[13]-1] == [7, 8, 11, 15, 16]
        @test vgraph.edge_data[off[13]:off[14]-1] == [9, 10, 14]
        @test vgraph.edge_data[off[14]:off[15]-1] == [9, 10, 11, 13, 15]
        @test vgraph.edge_data[off[15]:off[16]-1] == [10, 11, 12, 14, 16]
        @test vgraph.edge_data[off[16]:off[17]-1] == [11, 12, 15]
    end

    @testset "vertex connectivity graph for a 4×4 element quad mesh" begin
        mesh = equispaced_rectangular_mesh(limits..., 4, 4, (true, true))
        vgraph = build_vertex_graph(mesh)
        off = vgraph.edge_offset

        @test vgraph.nverts == 4 * 4
        @test vgraph.edge_data[off[1]:off[2]-1] == [2, 4, 5, 6, 8, 13, 14, 16]
        @test vgraph.edge_data[off[2]:off[3]-1] == [1, 3, 5, 6, 7, 13, 14, 15]
        @test vgraph.edge_data[off[3]:off[4]-1] == [2, 4, 6, 7, 8, 14, 15, 16]
        @test vgraph.edge_data[off[4]:off[5]-1] == [1, 3, 5, 7, 8, 13, 15, 16]
        @test vgraph.edge_data[off[5]:off[6]-1] == [1, 2, 4, 6, 8, 9, 10, 12]
        @test vgraph.edge_data[off[6]:off[7]-1] == [1, 2, 3, 5, 7, 9, 10, 11]
        @test vgraph.edge_data[off[7]:off[8]-1] == [2, 3, 4, 6, 8, 10, 11, 12]
        @test vgraph.edge_data[off[8]:off[9]-1] == [1, 3, 4, 5, 7, 9, 11, 12]
        @test vgraph.edge_data[off[9]:off[10]-1] == [5, 6, 8, 10, 12, 13, 14, 16]
        @test vgraph.edge_data[off[10]:off[11]-1] == [5, 6, 7, 9, 11, 13, 14, 15]
        @test vgraph.edge_data[off[11]:off[12]-1] == [6, 7, 8, 10, 12, 14, 15, 16]
        @test vgraph.edge_data[off[12]:off[13]-1] == [5, 7, 8, 9, 11, 13, 15, 16]
        @test vgraph.edge_data[off[13]:off[14]-1] == [1, 2, 4, 9, 10, 12, 14, 16]
        @test vgraph.edge_data[off[14]:off[15]-1] == [1, 2, 3, 9, 10, 11, 13, 15]
        @test vgraph.edge_data[off[15]:off[16]-1] == [2, 3, 4, 10, 11, 12, 14, 16]
        @test vgraph.edge_data[off[16]:off[17]-1] == [1, 3, 4, 9, 11, 12, 13, 15]
    end
end

@testset "graph for a cube panel mesh" begin
    @testset "face connectivity graph for a 6-element cube panel mesh" begin
        mesh = cube_panel_mesh(1, FT)
        fgraph = build_face_graph(mesh)
        off = fgraph.edge_offset

        @test fgraph.nverts == 6
        @test fgraph.edge_data[off[1]:off[2]-1] == [2, 3, 4, 5]
        @test fgraph.edge_data[off[2]:off[3]-1] == [1, 3, 5, 6]
        @test fgraph.edge_data[off[3]:off[4]-1] == [1, 2, 4, 6]
        @test fgraph.edge_data[off[4]:off[5]-1] == [1, 3, 5, 6]
        @test fgraph.edge_data[off[5]:off[6]-1] == [1, 2, 4, 6]
        @test fgraph.edge_data[off[6]:off[7]-1] == [2, 3, 4, 5]
    end

    @testset "vertex connectivity graph for a 6-element cube panel mesh" begin
        mesh = cube_panel_mesh(1, FT)
        vgraph = build_vertex_graph(mesh)
        off = vgraph.edge_offset

        @test vgraph.nverts == 6
        @test vgraph.edge_data[off[1]:off[2]-1] == [2, 3, 4, 5]
        @test vgraph.edge_data[off[2]:off[3]-1] == [1, 3, 5, 6]
        @test vgraph.edge_data[off[3]:off[4]-1] == [1, 2, 4, 6]
        @test vgraph.edge_data[off[4]:off[5]-1] == [1, 3, 5, 6]
        @test vgraph.edge_data[off[5]:off[6]-1] == [1, 2, 4, 6]
        @test vgraph.edge_data[off[6]:off[7]-1] == [2, 3, 4, 5]
    end

    @testset "face connectivity graph for 24 element cube panel mesh with 4 elements per panel" begin
        mesh = cube_panel_mesh(2, FT)
        fgraph = build_face_graph(mesh)
        off = fgraph.edge_offset

        @test fgraph.nverts == 24
        @test fgraph.edge_data[off[1]:off[2]-1] == [2, 3, 13, 17]
    end

    @testset "face connectivity graph for 24 element cube panel mesh with 4 elements per panel" begin
        mesh = cube_panel_mesh(2, FT)
        vgraph = build_vertex_graph(mesh)
        off = vgraph.edge_offset

        @test vgraph.nverts == 24
        @test vgraph.edge_data[off[1]:off[2]-1] == [2, 3, 4, 13, 15, 17, 18]
    end
end

@testset "graph laplacian matrix x vector" begin
    @testset "for a 4×4 element quad mesh" begin
        mesh = equispaced_rectangular_mesh(limits..., 4, 4, (false, false))
        fgraph = build_face_graph(mesh)
        nverts = fgraph.nverts
        vin = ones(FT, nverts)
        vout = Vector{FT}(undef, nverts)
        laplacian_x_v!(vout, vin, fgraph)
        @test vout == zeros(FT, nverts)
    end

    @testset "for a 24 element cube panel mesh" begin
        mesh = cube_panel_mesh(2, FT)
        fgraph = build_face_graph(mesh)
        nverts = fgraph.nverts
        vin = ones(FT, nverts)
        vout = Vector{FT}(undef, nverts)
        laplacian_x_v!(vout, vin, fgraph)
        @test vout == zeros(FT, nverts)
    end
end

@testset "test Fiedler vector calculation using iterative methods" begin
    @testset "for a 4×4 element quad mesh" begin
        mesh = equispaced_rectangular_mesh(limits..., 4, 4, (false, false))
        graph = build_face_graph(mesh)
        fdlr_direct  = fiedler_vector(graph, Graphs.Direct())
        fdlr_power   = fiedler_vector(graph, Graphs.PowerIt())
        fdlr_lanczos = fiedler_vector(graph, Graphs.Lanczos())
        @test is_in_fiedler_vector_space(fdlr_direct, graph)
        @test is_in_fiedler_vector_space(fdlr_power, graph)
        @test is_in_fiedler_vector_space(fdlr_lanczos, graph)
    end
end
