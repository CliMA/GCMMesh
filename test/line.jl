using Test
using Plots

using GCMMesh
using GCMMesh.Mesh

FT = Float64

limits = (FT(0), FT(1))

@testset "simple line mesh opposing face" begin

    @testset "assert correct element numbering" begin
        mesh = equispaced_line_mesh(limits..., 1, false)
        @test_throws AssertionError Mesh.opposing_face(mesh, 0, 1)
        @test_throws AssertionError Mesh.opposing_face(mesh, 2, 1)
    end

    @testset "1 element line mesh with all periodic boundries" begin
        mesh = equispaced_line_mesh(limits..., 1, true)
        @test Mesh.opposing_face(mesh, 1, 1) == (1, 2)
        @test Mesh.opposing_face(mesh, 1, 2) == (1, 1)
    end

    @testset "2 element line mesh with non-periodic boundaries" begin
        mesh = equispaced_line_mesh(limits..., 2, false)
        @test Mesh.opposing_face(mesh, 1, 1) == (0, 1)
        @test Mesh.opposing_face(mesh, 1, 2) == (2, 1)
        @test Mesh.opposing_face(mesh, 2, 1) == (1, 2)
        @test Mesh.opposing_face(mesh, 2, 2) == (0, 2)
    end

    @testset "2 element line mesh with periodic boundaries" begin
        mesh = equispaced_line_mesh(limits..., 2, true)
        @test Mesh.opposing_face(mesh, 1, 1) == (2, 2)
        @test Mesh.opposing_face(mesh, 1, 2) == (2, 1)
        @test Mesh.opposing_face(mesh, 2, 1) == (1, 2)
        @test Mesh.opposing_face(mesh, 2, 2) == (1, 1)
    end
end

@testset "simple line grid interior faces iterator" begin
    @testset "1×1 element line mesh with periodic boundary" begin
        mesh = equispaced_line_mesh(limits..., 1, true)
        @test length(Mesh.interior_faces(mesh)) == 1
        faces = Mesh.interior_faces(mesh)
        @test tuple(mesh.vertex_neighbors[faces[1], :]...) == (1, 2, 1, 1)
    end
    @testset "1×1 element quad mesh with non-periodic boundaries" begin
        mesh = equispaced_line_mesh(limits..., 1, false)
        @test Mesh.interior_faces(mesh) == nothing
    end
    @testset "2 element line mesh with non-periodic boundaries" begin
        mesh = equispaced_line_mesh(limits..., 2, false)
        @test length(Mesh.interior_faces(mesh)) == 1
        faces = Mesh.interior_faces(mesh)
        @test tuple(mesh.vertex_neighbors[faces[1], :]...) == (1, 2, 2, 1)
    end
end
@testset "simple line grid boundary faces iterator" begin
    @testset "1×1 element line mesh with all periodic boundries" begin
        mesh = equispaced_line_mesh(limits..., 1, true)
        @test Mesh.boundary_faces(mesh, 1) == nothing
        @test Mesh.boundary_faces(mesh, 2) == nothing
    end
    @testset "1×1 element quad mesh with non-periodic boundaries" begin
        mesh = equispaced_line_mesh(limits..., 1, false)
        faces1 = Mesh.boundary_faces(mesh, 1)
        faces2 = Mesh.boundary_faces(mesh, 2)

        @test length(faces1) == 1
        @test length(faces2) == 1
        @test [Mesh.boundary_face_info(mesh, face) for face in faces1] == [(1, 1)]
        @test [Mesh.boundary_face_info(mesh, face) for face in faces2] == [(1, 2)]
        @test Mesh.boundary_tag(mesh, :left) == 1
        @test Mesh.boundary_tag(mesh, :right) == 2
    end
    @testset "2 element line mesh with non-periodic boundaries" begin
        mesh = equispaced_line_mesh(limits..., 2, false)
        faces1 = Mesh.boundary_faces(mesh, 1)
        faces2 = Mesh.boundary_faces(mesh, 2)
        @test length(faces1) == 1
        @test length(faces2) == 1

        @test [Mesh.boundary_face_info(mesh, face) for face in faces1] == [(1, 1)]
        @test [Mesh.boundary_face_info(mesh, face) for face in faces2] == [(2, 2)]
    end
end

@testset "simple line grid vertex iterator" begin
    @testset "1×1 element line mesh with all periodic boundries" begin
        mesh = equispaced_line_mesh(limits..., 1, true)
        uverts = mesh.unique_verts
        @test length(uverts) == 1
        @test sort(Mesh.unique_vertex_connectivity(mesh, 1)) == sort([(1, 1), (1, 2)])
        @test length(Mesh.unique_vertex_connectivity(mesh, 1)) == 2
    end
    @testset "1×1 element line mesh with non-periodic boundaries" begin
        mesh = equispaced_line_mesh(limits..., 1, false)
        uverts = mesh.unique_verts
        @test length(uverts) == 2
        vc1 = Mesh.unique_vertex_connectivity(mesh, 1)
        vc2 = Mesh.unique_vertex_connectivity(mesh, 2)
        @test length(vc1) == 1
        @test collect(vc1) == [(1, 1)]
        @test length(vc2) == 1
        @test collect(vc2) == [(1, 2)]
    end
    @testset "2 element line mesh with non-periodic boundaries" begin
        mesh = equispaced_line_mesh(limits..., 2, false)
        uverts = mesh.unique_verts
        @test length(uverts) == 3
        vc1 = Mesh.unique_vertex_connectivity(mesh, 1)
        vc2 = Mesh.unique_vertex_connectivity(mesh, 2)
        vc3 = Mesh.unique_vertex_connectivity(mesh, 3)
        @test length(vc1) == 1
        @test vc1 == [(1, 1)]
        @test length(vc2) == 2
        @test vc2 == [(1, 2), (2, 1)]
        @test length(vc3) == 1
        @test vc3 == [(2, 2)]
    end
    @testset "2 element line mesh with periodic boundary" begin
        mesh = equispaced_line_mesh(limits..., 2, true)
        uverts = mesh.unique_verts
        @test length(uverts) == 2
        vc1 = Mesh.unique_vertex_connectivity(mesh, 1)
        vc2 = Mesh.unique_vertex_connectivity(mesh, 2)
        @test length(vc1) == 2
        @test sort(vc1) == [(1, 1), (2, 2)]
        @test length(vc2) == 2
        @test sort(vc2) == [(1, 2), (2, 1)]
    end
end

@testset "simple line grid coordinates" begin
    @testset "1 element line mesh with all periodic boundries" begin
        limits = (FT(0), FT(1))
        mesh = equispaced_line_mesh(limits..., 1, true)
        c1, c2 = Mesh.vertex_coordinates(mesh, 1)
        @test c1 == (0.0,)
        @test c2 == (1.0,)

        limits = (FT(-1), FT(1))
        mesh = equispaced_line_mesh(limits..., 1, false)
        c1, c2 = Mesh.vertex_coordinates(mesh, 1)
        @test c1 == (-1.0,)
        @test c2 == (1.0,)
    end

    @testset "2 element line mesh with non-periodic boundaries" begin
        limits = (FT(0), FT(1))
        mesh = equispaced_line_mesh(limits..., 2, false)
        c1, c2 = Mesh.vertex_coordinates(mesh, 1)
        @test c1 == (0.0,)
        @test c2 == (0.5,)

        c1, c2 = Mesh.vertex_coordinates(mesh, 2)
        @test c1 == (0.5,)
        @test c2 == (1.0,)
    end

    @testset "check coordinate type accuracy" begin
        limits = (big(0.0), big(1.0))
        mesh = equispaced_line_mesh(limits..., 3, false)
        c1, c2 = Mesh.vertex_coordinates(mesh, 1)
        @test eltype(c2) == BigFloat
        @test c2[1] ≈ big(1.0) / big(3.0) rtol = eps(BigFloat)
    end
end
