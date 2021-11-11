using Test
using Plots

using GCMMesh
using GCMMesh.Mesh

FT = Float64

x1min, x1max = FT(0), FT(1)
x2min, x2max = FT(0), FT(1)

limits = (x1min, x2min, x1max, x2max)

@testset "simple rectangular grid opposing face" begin

    @testset "assert correct element numbering" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (true, true))
        @test_throws AssertionError Mesh.opposing_face(mesh, 0, 1)
        @test_throws AssertionError Mesh.opposing_face(mesh, 2, 1)
    end

    @testset "1×1 element quad mesh with all periodic boundries" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (true, true))
        @test Mesh.opposing_face(mesh, 1, 1) == (1, 3, true)
        @test Mesh.opposing_face(mesh, 1, 2) == (1, 4, true)
        @test Mesh.opposing_face(mesh, 1, 3) == (1, 1, true)
        @test Mesh.opposing_face(mesh, 1, 4) == (1, 2, true)
    end

    @testset "1×1 element quad mesh with 1 periodic boundary" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (true, false))
        @test Mesh.opposing_face(mesh, 1, 1) == (0, 1, true)
        @test Mesh.opposing_face(mesh, 1, 2) == (1, 4, true)
        @test Mesh.opposing_face(mesh, 1, 3) == (0, 3, true)
        @test Mesh.opposing_face(mesh, 1, 4) == (1, 2, true)
    end

    @testset "1×1 element quad mesh with non-periodic boundaries" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (false, false))
        @test Mesh.opposing_face(mesh, 1, 1) == (0, 1, true)
        @test Mesh.opposing_face(mesh, 1, 2) == (0, 2, true)
        @test Mesh.opposing_face(mesh, 1, 3) == (0, 3, true)
        @test Mesh.opposing_face(mesh, 1, 4) == (0, 4, true)
    end

    @testset "2×2 element quad mesh with non-periodic boundaries" begin
        mesh = equispaced_rectangular_mesh(limits..., 2, 2, (false, false))
        @test Mesh.opposing_face(mesh, 1, 1) == (0, 1, true)
        @test Mesh.opposing_face(mesh, 1, 2) == (2, 4, true)
        @test Mesh.opposing_face(mesh, 1, 3) == (3, 1, true)
        @test Mesh.opposing_face(mesh, 1, 4) == (0, 4, true)
        @test Mesh.opposing_face(mesh, 2, 1) == (0, 1, true)
        @test Mesh.opposing_face(mesh, 2, 2) == (0, 2, true)
        @test Mesh.opposing_face(mesh, 2, 3) == (4, 1, true)
        @test Mesh.opposing_face(mesh, 2, 4) == (1, 2, true)
    end
end

@testset "simple rectangular grid interior faces iterator" begin
    @testset "1×1 element quad mesh with all periodic boundries" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (true, true))
        @test length(Mesh.interior_faces(mesh)) == 2
        faces = Mesh.interior_faces(mesh)
        @test tuple(mesh.face_neighbors[faces[1], :]...) == (1, 4, 1, 2, -1)
        @test tuple(mesh.face_neighbors[faces[2], :]...) == (1, 1, 1, 3, -1)
    end
    @testset "1×1 element quad mesh with 1 periodic boundary" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (true, false))
        @test length(Mesh.interior_faces(mesh)) == 1
        faces = Mesh.interior_faces(mesh)
        @test tuple(mesh.face_neighbors[faces[1], :]...) == (1, 4, 1, 2, -1)
    end
    @testset "1×1 element quad mesh with non-periodic boundaries" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (false, false))
        @test Mesh.interior_faces(mesh) == nothing
    end
    @testset "2×2 element quad mesh with non-periodic boundaries" begin
        mesh = equispaced_rectangular_mesh(limits..., 2, 2, (false, false))
        @test length(Mesh.interior_faces(mesh)) == 4
        faces = Mesh.interior_faces(mesh)
        @test tuple(mesh.face_neighbors[faces[1], :]...) == (1, 2, 2, 4, -1)
        @test tuple(mesh.face_neighbors[faces[2], :]...) == (3, 2, 4, 4, -1)
        @test tuple(mesh.face_neighbors[faces[3], :]...) == (1, 3, 3, 1, -1)
        @test tuple(mesh.face_neighbors[faces[4], :]...) == (2, 3, 4, 1, -1)
    end
end

@testset "simple rectangular grid boundary faces iterator" begin
    @testset "1×1 element quad mesh with all periodic boundries" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (true, true))
        @test Mesh.boundary_faces(mesh, 1) == nothing
        @test Mesh.boundary_faces(mesh, 2) == nothing
        @test Mesh.boundary_faces(mesh, 3) == nothing
        @test Mesh.boundary_faces(mesh, 4) == nothing
    end
    @testset "1×1 element quad mesh with 1 periodic boundary" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (true, false))
        @test Mesh.boundary_faces(mesh, 1) == nothing
        @test Mesh.boundary_faces(mesh, 2) == nothing
        faces3 = Mesh.boundary_faces(mesh, 3)
        faces4 = Mesh.boundary_faces(mesh, 4)
        @test length(faces3) == 1
        @test length(faces4) == 1
        @test [Mesh.boundary_face_info(mesh, face) for face in faces3] == [(1, 1)]
        @test [Mesh.boundary_face_info(mesh, face) for face in faces4] == [(1, 3)]
    end
    @testset "1×1 element quad mesh with non-periodic boundaries" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (false, false))
        faces1 = Mesh.boundary_faces(mesh, 1)
        faces2 = Mesh.boundary_faces(mesh, 2)
        faces3 = Mesh.boundary_faces(mesh, 3)
        faces4 = Mesh.boundary_faces(mesh, 4)

        @test length(faces1) == 1
        @test length(faces2) == 1
        @test length(faces3) == 1
        @test length(faces4) == 1
        @test [Mesh.boundary_face_info(mesh, face) for face in faces1] == [(1, 4)]
        @test [Mesh.boundary_face_info(mesh, face) for face in faces2] == [(1, 2)]
        @test [Mesh.boundary_face_info(mesh, face) for face in faces3] == [(1, 1)]
        @test [Mesh.boundary_face_info(mesh, face) for face in faces4] == [(1, 3)]
        @test Mesh.boundary_tag(mesh, :west) == 1
        @test Mesh.boundary_tag(mesh, :east) == 2
        @test Mesh.boundary_tag(mesh, :south) == 3
        @test Mesh.boundary_tag(mesh, :north) == 4
        @test Mesh.boundary_tag(mesh, :northeast) == nothing # no such boundary for rectangle
    end
    @testset "2×3 element quad mesh with non-periodic boundaries" begin
        mesh = equispaced_rectangular_mesh(limits..., 2, 3, (false, false))
        faces1 = Mesh.boundary_faces(mesh, 1)
        faces2 = Mesh.boundary_faces(mesh, 2)
        faces3 = Mesh.boundary_faces(mesh, 3)
        faces4 = Mesh.boundary_faces(mesh, 4)
        @test length(faces1) == 3
        @test length(faces2) == 3
        @test length(faces3) == 2
        @test length(faces4) == 2

        @test [Mesh.boundary_face_info(mesh, face) for face in faces1] == [(1, 4), (3, 4), (5, 4)]
        @test [Mesh.boundary_face_info(mesh, face) for face in faces2] == [(2, 2), (4, 2), (6, 2)]
        @test [Mesh.boundary_face_info(mesh, face) for face in faces3] == [(1, 1), (2, 1)]
        @test [Mesh.boundary_face_info(mesh, face) for face in faces4] == [(5, 3), (6, 3)]
    end
end


@testset "simple rectangular grid vertex iterator" begin
    @testset "1×1 element quad mesh with all periodic boundries" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (true, true))
        uverts = mesh.unique_verts
        @test length(uverts) == 1
        @test Mesh.unique_vertex_connectivity(mesh, 1) == [(1, 1), (1, 2), (1, 3), (1, 4)]
        @test length(Mesh.unique_vertex_connectivity(mesh, 1)) == 4
    end
    @testset "1×1 element quad mesh with 1 periodic boundary" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (true, false))
        uverts = mesh.unique_verts
        @test length(uverts) == 2
        vc1 = Mesh.unique_vertex_connectivity(mesh, 1)
        vc2 = Mesh.unique_vertex_connectivity(mesh, 2)
        @test length(vc1) == 2
        @test vc1 == [(1, 1), (1, 2)]
        @test length(vc2) == 2
        @test vc2 == [(1, 3), (1, 4)]
    end
    @testset "1×1 element quad mesh with non-periodic boundaries" begin
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (false, false))
        uverts = mesh.unique_verts
        @test length(uverts) == 4
        vc1 = Mesh.unique_vertex_connectivity(mesh, 1)
        vc2 = Mesh.unique_vertex_connectivity(mesh, 2)
        vc3 = Mesh.unique_vertex_connectivity(mesh, 3)
        vc4 = Mesh.unique_vertex_connectivity(mesh, 4)
        @test length(vc1) == 1
        @test collect(vc1) == [(1, 1)]
        @test length(vc2) == 1
        @test collect(vc2) == [(1, 2)]
        @test length(vc3) == 1
        @test collect(vc3) == [(1, 4)]
        @test length(vc4) == 1
        @test collect(vc4) == [(1, 3)]
    end
    @testset "2×3 element quad mesh with non-periodic boundaries" begin
        mesh = equispaced_rectangular_mesh(limits..., 2, 3, (false, false))
        uverts = mesh.unique_verts
        @test length(uverts) == 3 * 4
        vc1 = Mesh.unique_vertex_connectivity(mesh, 1)
        @test length(vc1) == 1
        @test vc1 == [(1, 1)]
    end
    @testset "2×3 element quad mesh with periodic boundaries" begin
        mesh = equispaced_rectangular_mesh(limits..., 2, 3, (true, true))
        uverts = mesh.unique_verts
        @test length(uverts) == 2 * 3
        vc1 = Mesh.unique_vertex_connectivity(mesh, 1)
        vc2 = Mesh.unique_vertex_connectivity(mesh, 2)
        vc3 = Mesh.unique_vertex_connectivity(mesh, 3)
        vc4 = Mesh.unique_vertex_connectivity(mesh, 4)
        vc5 = Mesh.unique_vertex_connectivity(mesh, 5)
        vc6 = Mesh.unique_vertex_connectivity(mesh, 6)
        @test sort(vc1) == sort([(1, 1), (2, 2), (5, 4), (6, 3)])
        @test sort(vc2) == sort([(2, 1), (1, 2), (5, 3), (6, 4)])
        @test sort(vc3) == sort([(3, 1), (4, 2), (1, 4), (2, 3)])
        @test sort(vc4) == sort([(4, 1), (3, 2), (2, 4), (1, 3)])
        @test sort(vc5) == sort([(5, 1), (6, 2), (3, 4), (4, 3)])
        @test sort(vc6) == sort([(6, 1), (5, 2), (4, 4), (3, 3)])
    end
end

@testset "simple rectangular grid coordinates" begin
    @testset "1×1 element quad mesh with all periodic boundries" begin
        x1min, x1max = FT(0), FT(1)
        x2min, x2max = FT(0), FT(1)
        limits = (x1min, x2min, x1max, x2max)
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (true, true))
        c1, c2, c3, c4 = Mesh.vertex_coordinates(mesh, 1)
        @test c1 == (0.0, 0.0)
        @test c2 == (1.0, 0.0)
        @test c3 == (1.0, 1.0)
        @test c4 == (0.0, 1.0)

        x1min, x1max = FT(-1), FT(1)
        x2min, x2max = FT(-1), FT(1)
        limits = (x1min, x2min, x1max, x2max)
        mesh = equispaced_rectangular_mesh(limits..., 1, 1, (false, false))
        c1, c2, c3, c4 = Mesh.vertex_coordinates(mesh, 1)
        @test c1 == (-1.0, -1.0)
        @test c2 == (1.0, -1.0)
        @test c3 == (1.0, 1.0)
        @test c4 == (-1.0, 1.0)
    end
    @testset "2×4 element quad mesh with non-periodic boundaries" begin
        x1min, x1max = FT(0), FT(1)
        x2min, x2max = FT(0), FT(1)
        limits = (x1min, x2min, x1max, x2max)
        mesh = equispaced_rectangular_mesh(limits..., 2, 4, (false, false))
        c1, c2, c3, c4 = Mesh.vertex_coordinates(mesh, 1)
        @test c1 == (0.0, 0.0)
        @test c2 == (0.5, 0.0)
        @test c3 == (0.5, 0.25)
        @test c4 == (0.0, 0.25)

        c1, c2, c3, c4 = Mesh.vertex_coordinates(mesh, 8)
        @test c1 == (0.5, 0.75)
        @test c2 == (1.0, 0.75)
        @test c3 == (1.0, 1.0)
        @test c4 == (0.5, 1.0)
    end
    @testset "check coordinate type accuracy" begin
        x1min, x1max = big(0.0), big(1.0)
        x2min, x2max = big(0.0), big(1.0)
        limits = (x1min, x2min, x1max, x2max)
        mesh = equispaced_rectangular_mesh(limits..., 3, 1, (false, false))
        c1, c2, c3, c4 = Mesh.vertex_coordinates(mesh, 1)
        @test eltype(c2) == BigFloat
        @test c2[1] ≈ big(1.0) / big(3.0) rtol = eps(BigFloat)
        @test c2[2] == 0.0
    end
end
