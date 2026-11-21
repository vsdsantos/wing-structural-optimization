"""
Regression tests for the structural-analysis path.

The original MATLAB line-integral helpers apply `trapz` to segment-weighted
values with implicit unit spacing. Julia instead uses the segment midpoint
sum, which is well-defined for a single segment and exact for the supported
piecewise-constant wall model. Analytical checks use tight floating-point
tolerances; load integration uses 1e-10 relative tolerance. For simply
supported loads, Julia resolves reactions at the nearest sampled supports and
uses root-referenced equilibrium; off-grid support results may therefore
differ from MATLAB's mixed-coordinate calculation.
"""

using Test
using WingOpt

function test_section()
    inertia = [2.0 0.0; 0.0 1.0]
    return SectionProperties(
        inertia, inertia, 2.0, 1.0, 0.0, [0.2, 0.0], 0.1, 0.02,
        0.5, 0.2, 0.18, [-0.05 0.1; 0.45 0.1; -0.05 -0.1; 0.45 -0.1],
        [-0.045 0.095; 0.445 0.095; -0.045 -0.095; 0.445 -0.095],
        [0.2, 0.095], [0.2, -0.095],
    )
end

@testset "Structural numerical primitives" begin
    @test trapz([0.0, 1.0], [1.0, 3.0]) ≈ 2.0
    @test trapz(Float64[]) == 0.0
    @test cumtrapz([0.0, 1.0, 2.0], [0.0, 1.0, 2.0]) ≈ [0.0, 0.5, 2.0]
    @test_throws DimensionMismatch trapz([0.0], [1.0, 2.0])
    @test sort(collect(quadratic_roots(1.0, -3.0, 2.0))) ≈ [1.0, 2.0]

    translated = translate_inertia_tensor([1.0 0.0; 0.0 2.0], [2.0, 3.0], 4.0)
    @test translated ≈ [17.0 24.0; 24.0 38.0]
    centroid = section_centroid([0.0 0.0; 2.0 0.0])
    @test centroid ≈ [1.0, 0.0]
    @test area_inertia_tensor(0.5, [0.0 2.0; 2.0 2.0]) ≈ [4.0 2.0; 2.0 1.0]
    @test multielement_centroid([0.0 0.0; 2.0 0.0], [1.0, 3.0]) ≈ [1.5, 0.0]
end

@testset "Section, material, and failure criteria" begin
    airfoil = parse_foil("s1223"; foil_dir = joinpath(@__DIR__, "../../foils"))
    foil_geometry = foil_properties(1.0, 0.0, [0.2, 0.7], [0.0, 0.0], airfoil)
    @test foil_geometry.area > 0
    @test foil_geometry.perimeter > 0

    foil = FoilData(x -> 0.1, x -> -0.1, zeros(2, 2), zeros(2, 2), 0.25)
    station = Station(1.0, 0.01, 0.005, 0.0, [0.25, 0.75])
    props = reinforced_box_section(station, foil)
    @test props.width ≈ 0.5
    @test props.height ≈ 0.2
    @test props.area ≈ 0.0136
    @test props.area_mid ≈ 0.0931
    @test props.Iz0 > props.Iy0 > 0

    balsa = material_properties("balsa")
    @test balsa.dens == 150.0
    @test_throws ArgumentError material_properties("unknown")
    stress = zeros(3, 3)
    stress[1, 1] = 20e6
    @test max_tension_criteria(stress, balsa, "axial") ≈ 1.0
    @test tsai_hill_criteria(stress, balsa) ≈ 0.25
    @test tsai_wu_criteria(stress, balsa, "axial") > 0
    @test isinf(max_tension_criteria(zeros(3, 3), balsa, "axial"))
    station = Station(1.0, 0.01, 0.005, 0.0, [0.25, 0.75])
    loads = HalfLoads([0.0], [100.0], [10.0], [1.0], [0.0], 1)
    @test isfinite(safety_margin(station, loads, 1, test_section()))
end

@testset "Thin-wall torsion and beam bending" begin
    phi, dphi = ang_torcao_parede_fina(5.0, [0.01, 0.01], [1.0, 1.0], 0.1, 10.0, [0.0, 2.0])
    @test integral_linha_media([0.01, 0.01], [1.0, 1.0]) ≈ 200.0
    @test dphi(0.5) ≈ 2500.0
    @test phi ≈ 5000.0
    @test_throws DomainError integral_linha_media([0.0], [1.0])

    props = test_section()
    stations = [Station(1.0, 0.01, 0.005, s, [0.25, 0.75]) for s in 0.0:1.0:3.0]
    half = HalfLoads(
        collect(0.0:1.0:3.0), zeros(4), fill(10.0, 4), fill(5.0, 4), zeros(4), 1
    )
    loads = InternalLoads(half, half)
    mat = MaterialProperties(1e8, 1e8, 100.0, 0.3, 1e8, 1e8, 100.0, 0.3, 10.0, 1e8, 100.0)
    props_vec = fill(props, 4)
    @test wing_torsion(stations, props_vec, [0.0, 5.0, 5.0, 5.0], mat) ≈ [
        0.0, 875.0, 2625.0, 4375.0
    ]
    @test wing_displacement(stations, props_vec, loads, mat).R ≈ [0.0, 0.05, 0.1, 0.1]
    report_wing = WingGeometry(
        1.0, 6.0, [1.0, 0.5], 0.75, [0.0, 0.0], [0.0, 0.0], [0.0, 0.0],
        ["test", "test"], 1.0, 0.0, [1.0, 0.5], [0.0, 1.0, 3.0],
    )
    @test isfinite(report_45(stations, props_vec, report_wing, mat; design_speed = 10.0))
    irregular = [
        stations[1], stations[2], Station(1.0, 0.01, 0.005, 2.5, [0.25, 0.75]),
        stations[4],
    ]
    @test_throws ArgumentError wing_displacement(irregular, props_vec, loads, mat)
end

@testset "Simply-supported wing internal loads and mass" begin
    b = 4.0
    wing = WingGeometry(
        1.0, b, [1.0, 0.5], 0.75, [0.0, 0.0], [0.0, 0.0], [0.0, 0.0],
        ["test", "test"], 1.0, 0.0, [1.0, 0.5], [0.0, b / 2],
    )
    aero = AeroLoad([-1.5, -0.5, 0.5, 1.5], [1.0, 2.0, 2.0, 1.0], zeros(4))
    props = fill(test_section(), 2)
    internal = internal_loads_SS(aero, props, wing)
    @test length(internal.L.s) == length(internal.R.s) == 2
    @test internal.L.s ≈ internal.R.s
    @test all(isfinite, internal.R.V)
    @test internal.L.V ≈ internal.R.V
    @test internal.L.M ≈ internal.R.M
    @test internal.L.Mt ≈ internal.R.Mt
    @test compute_internal_loads(aero, props, wing, "ss").R.M ≈ internal.R.M
    @test_throws ArgumentError compute_internal_loads(aero, props, wing, "invalid")

    stations = [Station(1.0, 0.01, 0.005, 0.0, [0.25, 0.75]),
                Station(1.0, 0.01, 0.005, 2.0, [0.25, 0.75])]
    mass, dm = beam_mass(
        stations, fill(test_section(), 2), material_properties("balsa"),
        material_properties("freijo"),
    )
    @test dm ≈ [4.2362]
    @test mass ≈ 4.2362
end
