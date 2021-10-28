using Test

using GCMMesh
using GCMMesh.Mesh

FT = Float64

@testset "simple cube surface grid opposing face" begin

    @testset "assert correct element numbering" begin
        mesh = cube_panel_mesh(1, FT)
        @test_throws AssertionError Mesh.opposing_face(
            mesh,
            0,
            1,
        )
        @test_throws AssertionError Mesh.opposing_face(
            mesh,
            7,
            1,
        )
    end
    @testset "6 element mesh with 1 element per panel" begin
        mesh = cube_panel_mesh(1, FT)
        @test Mesh.opposing_face(mesh, 1, 1) == (4, 4, true)
        @test Mesh.opposing_face(mesh, 1, 2) == (3, 1, true)
        @test Mesh.opposing_face(mesh, 1, 3) == (2, 1, true)
        @test Mesh.opposing_face(mesh, 1, 4) == (5, 1, true)

        @test Mesh.opposing_face(mesh, 2, 1) == (1, 3, true)
        @test Mesh.opposing_face(mesh, 2, 2) == (3, 4, true)
        @test Mesh.opposing_face(mesh, 2, 3) == (6, 2, true)
        @test Mesh.opposing_face(mesh, 2, 4) == (5, 2, true)

        @test Mesh.opposing_face(mesh, 3, 1) == (1, 2, true)
        @test Mesh.opposing_face(mesh, 3, 2) == (4, 3, true)
        @test Mesh.opposing_face(mesh, 3, 3) == (6, 3, true)
        @test Mesh.opposing_face(mesh, 3, 4) == (2, 2, true)

        @test Mesh.opposing_face(mesh, 4, 1) == (5, 4, true)
        @test Mesh.opposing_face(mesh, 4, 2) == (6, 4, true)
        @test Mesh.opposing_face(mesh, 4, 3) == (3, 2, true)
        @test Mesh.opposing_face(mesh, 4, 4) == (1, 1, true)

        @test Mesh.opposing_face(mesh, 5, 1) == (1, 4, true)
        @test Mesh.opposing_face(mesh, 5, 2) == (2, 4, true)
        @test Mesh.opposing_face(mesh, 5, 3) == (6, 1, true)
        @test Mesh.opposing_face(mesh, 5, 4) == (4, 1, true)

        @test Mesh.opposing_face(mesh, 6, 1) == (5, 3, true)
        @test Mesh.opposing_face(mesh, 6, 2) == (2, 3, true)
        @test Mesh.opposing_face(mesh, 6, 3) == (3, 3, true)
        @test Mesh.opposing_face(mesh, 6, 4) == (4, 2, true)

        # 6 faces, 4 vertices per face, 8 global vertices, so each vertex should be part of 3 elements        
        # check that all vertices appear as part of 3 elements
        for uvertno in 1:length(mesh.unique_verts)
            @test length(Mesh.unique_vertex_connectivity(mesh, uvertno)) == 3
        end
    end
    @testset "24 element mesh with 4 elements per panel" begin
        mesh = cube_panel_mesh(2, FT)
        @test Mesh.opposing_face(mesh, 1, 1) == (13, 4, true)
        @test Mesh.opposing_face(mesh, 1, 2) == (2, 4, true)
        @test Mesh.opposing_face(mesh, 1, 3) == (3, 1, true)
        @test Mesh.opposing_face(mesh, 1, 4) == (17, 1, true)
    end
end

@testset "cube surface grid interior faces iterator" begin
    @testset "all faces should be interior faces" begin
        mesh = cube_panel_mesh(1, FT)
        @test length(Mesh.interior_faces(mesh)) ==
              mesh.nfaces
        mesh = cube_panel_mesh(2, FT)
        @test length(Mesh.interior_faces(mesh)) ==
              mesh.nfaces
    end

    @testset "6 element mesh with 1 element per panel" begin
        mesh = cube_panel_mesh(1, FT)
        faces = Mesh.interior_faces(mesh)
        @test tuple(mesh.face_neighbors[faces[1], :]...) == (1, 4, 5, 1, -1)
        @test tuple(mesh.face_neighbors[faces[2], :]...) == (1, 3, 2, 1, -1)
        @test tuple(mesh.face_neighbors[faces[3], :]...) == (1, 2, 3, 1, -1)
        @test tuple(mesh.face_neighbors[faces[4], :]...) == (1, 1, 4, 4, -1)
        @test tuple(mesh.face_neighbors[faces[5], :]...) == (4, 1, 5, 4, -1)
        @test tuple(mesh.face_neighbors[faces[6], :]...) == (5, 2, 2, 4, -1)
        @test tuple(mesh.face_neighbors[faces[7], :]...) == (2, 2, 3, 4, -1)
        @test tuple(mesh.face_neighbors[faces[8], :]...) == (4, 3, 3, 2, -1)
        @test tuple(mesh.face_neighbors[faces[9], :]...) == (5, 3, 6, 1, -1)
        @test tuple(mesh.face_neighbors[faces[10], :]...) == (6, 2, 2, 3, -1)
        @test tuple(mesh.face_neighbors[faces[11], :]...) == (6, 3, 3, 3, -1)
        @test tuple(mesh.face_neighbors[faces[12], :]...) == (4, 2, 6, 4, -1)
    end
end

@testset "simple cube surface grid boundary faces iterator" begin
    @testset "24 element mesh with 4 elements per panel (no boundary faces for this topology)" begin
        mesh = cube_panel_mesh(2, FT)
        @test Mesh.boundary_faces(mesh, 1) == nothing
        @test Mesh.boundary_faces(mesh, 2) == nothing
        @test Mesh.boundary_faces(mesh, 3) == nothing
        @test Mesh.boundary_faces(mesh, 4) == nothing
    end
end

@testset "sphere mesh" begin
    FT = Float64
    @testset "4 elements per edge, equidistant spherical mesh of radius 10; Float type = Float64" begin
        radius = FT(10)
        mesh = sphere_mesh(4, radius,  EquidistantSphereWarp())
        crad = abs.(sqrt.(sum(mesh.coordinates .^ 2, dims = 2)) .- radius)
        @test maximum(crad) ≤ 100 * eps(FT)
    end
    @testset "4 elements per edge, equiangular spherical mesh of radius 10; Float type = Float64" begin
        radius = FT(10)
        mesh = sphere_mesh(4, radius, EquiangularSphereWarp())
        crad = abs.(sqrt.(sum(mesh.coordinates .^ 2, dims = 2)) .- radius)
        @test maximum(crad) ≤ 100 * eps(FT)
    end

    FT = BigFloat
    @testset "4 elements per edge, equidistant spherical mesh of radius 10; Float type = BigFloat" begin
        radius = FT(10)
        mesh = sphere_mesh(4, radius, EquidistantSphereWarp())
        crad = abs.(sqrt.(sum(mesh.coordinates .^ 2, dims = 2)) .- radius)
        @test maximum(crad) ≤ 100 * eps(FT)
    end
    @testset "4 elements per edge, equiangular spherical mesh of radius 10; Float type = BigFloat" begin
        radius = FT(10)
        mesh = sphere_mesh(4, radius, EquiangularSphereWarp())
        crad = abs.(sqrt.(sum(mesh.coordinates .^ 2, dims = 2)) .- radius)
        @test maximum(crad) ≤ 100 * eps(FT)
    end
end
