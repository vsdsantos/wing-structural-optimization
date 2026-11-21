"""
Unit tests for WingOpt.jl
Run with:  julia --project=WingOpt WingOpt/test/runtests.jl
or:        cd WingOpt && julia --project -e 'using Pkg; Pkg.test()'
"""

using Test
using LinearAlgebra
using WingOpt

@testset "Core numerical foundations" begin
    @testset "numerical integration and roots" begin
        x = range(0.0, 1.0; length = 101)
        y = x .^ 2
        @test trapz(x, y) ≈ 1 / 3 atol = 5e-5
        @test trapz([1.0, 2.0, 3.0]) == 4.0
        @test cumtrapz([0.0, 1.0, 2.0], [1.0, 1.0, 1.0]) == [0.0, 1.0, 2.0]
        @test cumtrapz(Float64[]) == Float64[]
        @test_throws DimensionMismatch trapz([0.0], [1.0, 2.0])
        @test_throws ArgumentError cumtrapz([0.0, Inf], [1.0, 2.0])

        roots = quadratic_roots(1.0, -3.0, 2.0)
        @test collect(roots) ≈ [2.0 + 0im, 1.0 + 0im]
        complex_roots = quadratic_roots(1.0, 0.0, 1.0)
        @test collect(complex_roots) ≈ [-1im, 1im]
        @test_throws DomainError quadratic_roots(0.0, 1.0, 1.0)
    end

    @testset "numerical properties" begin
        for n = 2:12
            x = range(-1.0, 2.0; length = n)
            y = sin.(x) .+ 2
            @test cumtrapz(x, y)[end] ≈ trapz(x, y)
            @test trapz(reverse(x), reverse(y)) ≈ -trapz(x, y)
        end
        inertia = [4.0 0.75; 0.75 2.0]
        for angle in range(-pi, pi; length = 9)
            rotated = rotate_inertia_tensor(inertia, angle)
            @test tr(rotated) ≈ tr(inertia)
            @test det(rotated) ≈ det(inertia)
        end
    end

    @testset "domain types and materials" begin
        wing = WingGeometry(
            1.0,
            7.0,
            [1.0, 0.5],
            0.75,
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            ["s1223", "s1223"],
            0.1,
            0.0,
            [1.0, 1.0, 0.5],
            [0.0, 2.625, 3.5],
        )
        @test wing.full_span == 7.0
        @test material_properties(" Balsa ").dens == 150
        @test material_properties("carbono_uni").axial_trac == 1.13e9
        @test_throws ArgumentError material_properties("unknown")
        @test_throws DomainError FoilGeometry(
            [0 0; 1 0],
            [0 0; 1 0],
            1,
            1,
            0,
            0,
            2,
            0,
            [0, 1],
            [0, 0],
        )
        @test_throws DomainError MaterialProperties(1, 1, 1, 0.5, 1, 1, 1, 0.2, 1, 1, 1)
        @test_throws DimensionMismatch WingGeometry(
            1,
            7,
            [1, 0.5],
            0.75,
            [0, 0],
            [0, 0],
            [0, 0],
            ["foil"],
            0.1,
            0,
            [1, 1, 0.5],
            [0, 2.625, 3.5],
        )
        @test_throws DomainError Station(1.0, 0.001, 0.0, 1.0, [0.8, 0.2])
        @test AeroLoad([0.0, 1.0], [1.0, 1.0], [0.0, 0.0]).dLdb == [1.0, 1.0]
        @test_throws DimensionMismatch AeroLoad([0.0], [1.0, 2.0], [0.0])
        @test_throws ArgumentError AeroLoad([1.0, 0.0], [1.0, 1.0], [0.0, 0.0])
    end

    @testset "section geometry and inertia transforms" begin
        points = [0.0 0.0; 1.0 0.0]
        @test section_centroid(points) == [0.5, 0.0]
        @test area_inertia_tensor(0.2, points) ≈ [0.0 0.0; 0.0 0.2 / 3]
        @test first_moment_of_area(points, 0.2) == (0.0, 0.1)
        @test multielement_centroid([0.0 0.0; 2.0 0.0], [1.0, 3.0]) == [1.5, 0.0]
        @test_throws DomainError multielement_centroid([0.0 0.0], [0.0])

        inertia = [3.0 0.5; 0.5 1.0]
        moved = translate_inertia_tensor(inertia, [2.0, -1.0], 4.0)
        @test moved == [19.0 -7.5; -7.5 5.0]
        principal, angle = principal_inertia_tensor(inertia)
        @test principal[1, 1] >= principal[2, 2]
        @test rotate_inertia_tensor(inertia, angle) ≈ principal
        @test principal[1, 1] + principal[2, 2] ≈ tr(inertia)
        @test_throws DimensionMismatch rotate_inertia_tensor(ones(3, 3), 0.0)
        @test_throws DomainError principal_inertia_tensor([-1.0 0.0; 0.0 1.0])
    end

    @testset "airfoil interpolation and spanwise stations" begin
        foil = parse_foil("s1223")
        @test foil.ext[1, 1] == minimum(foil.ext[:, 1])
        @test 0 <= foil.ext[1, 1] < foil.ext[end, 1] <= 1
        @test foil.ext[end, 1] == 1.0
        @test foil.esp > 0
        profile = foil_properties(1.0, 0.0, [0.1, 0.6], [0.0, 0.0], foil)
        @test size(profile.upper, 2) == 2
        @test profile.area > 0
        @test profile.perimeter > 0
        full_profile = foil_properties(1.0, 0.0, [0.0, 1.0], [0.0, 0.0], foil)
        @test full_profile.area > 0
        @test_throws DomainError foil_properties(1.0, 0.0, [-0.1, 0.6], [0.0, 0.0], foil)

        wing = WingGeometry(
            1.0,
            7.0,
            [1.0, 0.5],
            0.75,
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            ["s1223", "s1223"],
            0.1,
            0.0,
            [1.0, 1.0, 0.5],
            [0.0, 2.625, 3.5],
        )
        geometry = StructuredGeometry(
            [-0.15, -0.15, -0.1],
            [0.3, 0.3, 0.2],
            [0.002, 0.002, 0.001],
            [0.003, 0.003, 0.002],
        )
        stations = generate_wing_stations(wing, geometry, [0.0, 1.75, 3.5])
        @test length(stations) == 3
        @test stations[1].limits ≈ [0.1, 0.4]
        @test stations[end].chord == 0.5
        @test_throws DomainError generate_wing_stations(wing, geometry, [4.0])
    end
end
