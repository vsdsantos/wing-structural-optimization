"""
Unit tests for WingOpt.jl
Run with:  julia --project=WingOpt WingOpt/test/runtests.jl
or:        cd WingOpt && julia --project -e 'using Pkg; Pkg.test()'
"""

using Test
using WingOpt

const FOIL_DIR = normpath(joinpath(@__DIR__, "..", "..", "foils"))

@testset "Numerical foundations" begin
    @test WingOpt.trapz([0.0, 1.0], [0.0, 1.0]) ≈ 0.5
    @test WingOpt.cumtrapz([0.0, 1.0], [0.0, 1.0]) ≈ [0.0, 0.5]
    @test_throws DimensionMismatch WingOpt.trapz([0.0], [0.0, 1.0])
    @test collect(WingOpt.quadratic_roots(1.0, -3.0, 2.0)) ≈ [2.0, 1.0]

    translated = translate_inertia_tensor([1.0 0.0; 0.0 1.0], [2.0, 3.0], 4.0)
    @test translated ≈ [17.0 24.0; 24.0 37.0]
    principal, angle = principal_inertia_tensor([4.0 0.0; 0.0 2.0])
    @test principal ≈ [4.0 0.0; 0.0 2.0]
    @test angle ≈ 0.0
    @test material_properties("balsa").dens == 150.0
    @test_throws ArgumentError material_properties("not-a-material")
end

@testset "Configuration, airfoil, and section geometry" begin
    config = build_run_config(1)
    wing = build_wing_geometry(config)
    @test wing.span_transitions == [0.0, 2.625, 3.5]
    @test_throws ArgumentError build_run_config(2)

    foil = parse_foil(config.foil_name; foil_dir = FOIL_DIR)
    @test foil.ext[1, 1] ≈ 0.0
    @test foil.ext[end, 1] ≈ 1.0
    design = [0.2, 0.2, 0.2, 2.0, 0.005]
    geometry = reinforced_box_parametric_geometry(design, wing)
    stations = generate_wing_stations(wing, geometry, [0.5, 1.5, 3.0])
    @test length(stations) == 3
    @test all(station -> station.limits[1] < station.limits[2], stations)
    props = [reinforced_box_section(station, foil) for station in stations]
    @test all(prop -> prop.area > 0 && prop.Iz0 > 0 && prop.area_mid > 0, props)
    @test minimum(props[2].vertices[:, 1]) <=
          props[2].centroid[1] <=
          maximum(props[2].vertices[:, 1])
end

@testset "Failure criteria and optimization" begin
    material = material_properties("carbono_bi")
    zero_stress = zeros(3, 3)
    @test max_tension_criteria(zero_stress, material, "axial") == Inf
    @test tsai_wu_criteria(zero_stress, material, "axial") == Inf
    @test tsai_hill_criteria(zero_stress, material) == 1
    stress = zeros(3, 3)
    stress[1, 1] = material.axial_trac / 2
    @test max_tension_criteria(stress, material, "axial") > 0
    @test tsai_wu_criteria(stress, material, "axial") > 0
    @test weight_function(2.0, 2.0, 0.0) == -0.5

    result = run_optimizer(
        :ps,
        x -> sum(abs2, x .- 0.25);
        x0 = [0.2, 0.2, 0.2, 2.0, 0.005],
        max_iterations = 2,
    )
    @test result.fval isa Float64
    @test result.x[4] in 1.0:4.0
    @test result.x[1] + result.x[2] <= 0.75
end

@testset "Structural fitness smoke test" begin
    config = build_run_config(1)
    wing = build_wing_geometry(config)
    foil = parse_foil(config.foil_name; foil_dir = FOIL_DIR)
    loads = example_aero_loads(wing; n_panels = 8)
    fitness, data = fit_function([0.2, 0.2, 0.2, 2.0, 0.005], config, wing, foil, loads)
    @test isfinite(fitness)
    @test data.mass > 0
    @test length(data.stations) == 4
    @test length(data.internal_loads) == 1
    @test length(data.internal_loads[1].R.M) == length(data.stations)
    result, optimized_data = run_optimization(
        1;
        optimizer = :ps,
        max_iterations = 1,
        n_panels = 8,
        foil_dir = FOIL_DIR,
    )
    @test isfinite(result.fval)
    @test optimized_data.mass > 0
end
