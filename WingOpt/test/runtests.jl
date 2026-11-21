"""
Unit tests for WingOpt.jl
Run with:  julia --project=WingOpt WingOpt/test/runtests.jl
or:        cd WingOpt && julia --project -e 'using Pkg; Pkg.test()'
"""

using Test
using WingOpt

function test_wing(; span = 4.0, constraint_distance = 1.0, constraint_angle = 0.0)
    WingGeometry(
        1.0,
        span,
        [1.0, 0.5],
        0.75,
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0],
        ["test", "test"],
        constraint_distance,
        constraint_angle,
        [1.0, 0.5],
        [0.0, span / 2],
    )
end

function test_property(; centroid = [0.25, 0.0], area_mid = 0.5, width = 0.5, height = 0.5)
    SectionProperties(
        zeros(2, 2),
        zeros(2, 2),
        1.0,
        1.0,
        0.0,
        centroid,
        area_mid,
        1.0,
        width,
        height,
        0.5,
        zeros(4, 2),
        zeros(4, 2),
        zeros(2),
        zeros(2),
    )
end

function test_material(; G = 10.0)
    MaterialProperties(1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, G, 1.0, 1.0)
end

@testset "WingOpt" begin
    @testset "aerodynamic loads" begin
        wing = test_wing()
        loads = only(
            example_aero_loads(
                wing;
                n_panels = 10,
                total_lift = 20.0,
                total_pitching_moment = 4.0,
            ),
        )
        @test length(loads.x) == 10
        @test loads.x ≈ -reverse(loads.x)
        @test loads.dLdb ≈ reverse(loads.dLdb)
        @test WingOpt.trapz(vcat(-2.0, loads.x, 2.0), vcat(0.0, loads.dLdb, 0.0)) ≈ 20.0
        @test WingOpt.trapz(vcat(-2.0, loads.x, 2.0), vcat(0.0, loads.dMadb, 0.0)) ≈ 4.0
        @test_throws ArgumentError example_aero_loads(wing; n_panels = 1)
    end

    @testset "simply-supported internal loads" begin
        wing = test_wing(constraint_distance = 0.5)
        x = collect(range(-1.75; step = 0.5, length = 8))
        dLdb = ones(8)
        props = [
            test_property(centroid = [0.25 * (1 - 0.25 * s), 0.0]) for
            s in (0.25, 0.75, 1.25, 1.75)
        ]
        load = AeroLoad(x, dLdb, zeros(8))
        internal = compute_internal_loads(load, props, wing, "ss")
        @test length(internal.L.s) == length(internal.R.s) == 4
        @test issorted(internal.L.s) && issorted(internal.R.s)
        @test internal.L.V ≈ internal.R.V
        @test internal.L.M ≈ internal.R.M
        @test all(iszero, internal.R.F)
        @test all(iszero, internal.L.Mt) && all(iszero, internal.R.Mt)
        twisted = internal_loads_SS(AeroLoad(x, dLdb, ones(8)), props, wing)
        @test any(x -> !iszero(x), twisted.R.Mt)
        @test twisted.L.Mt ≈ twisted.R.Mt
        @test_throws ArgumentError compute_internal_loads(load, props, wing, "unknown")
        @test_throws DimensionMismatch internal_loads_SS(
            AeroLoad(x, dLdb[1:7], zeros(7)),
            props,
            wing,
        )
    end

    @testset "thin-wall torsion" begin
        line_integral = integral_linha_media([0.002, 0.003], [1.0, 2.0])
        @test line_integral ≈ 0.008
        phi, dphi_dx =
            ang_torcao_parede_fina(2.0, [0.002, 0.003], [1.0, 2.0], 1.0, 10.0, [0.0, 3.0])
        @test phi ≈ 0.0012
        @test dphi_dx(1.0) ≈ 0.0004
        @test integral_linha_media([0.002], [0.0]) == 0.0
        @test first(
            ang_torcao_parede_fina(
                x -> x,
                [0.002, 0.003],
                [1.0, 2.0],
                1.0,
                10.0,
                [0.0, 1.0],
            ),
        ) ≈ 0.0001
        @test first(ang_torcao_parede_fina(0.0, [0.002], [1.0], 1.0, 10.0, [2.0, 2.0])) ==
              0.0
        @test_throws DimensionMismatch integral_linha_media([1.0], [1.0, 2.0])
        @test_throws ArgumentError integral_linha_media([0.0], [1.0])
    end

    @testset "bi-cell and wing torsion" begin
        t_cells = [[1.0], [1.0], [1.0]]
        s_cells = [[1.0], [1.0], [1.0]]
        flux = bicelular([1.0, 2.0], [1.0, 1.0], [1.0, -1.0, 0.0], t_cells, s_cells)
        @test flux ≈ [0.75, 0.75, 0.0]
        stations = [Station(1.0, 0.01, 0.0, span, [0.1, 0.9]) for span in (0.0, 1.0)]
        props = [test_property() for _ in stations]
        @test wing_torsion(stations, props, [1.0, 1.0], test_material()) ≈ [0.0, 20.0]
        @test_throws DimensionMismatch wing_torsion(stations, props, [1.0], test_material())
    end
end
