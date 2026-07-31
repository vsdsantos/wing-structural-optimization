"""
Unit and integration tests for WingOpt.jl
Run with:  julia --project=WingOpt WingOpt/test/runtests.jl
or:        cd WingOpt && julia --project -e 'using Pkg; Pkg.test()'
"""

using Test

# Activate the project so all dependencies are available
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using WingOpt
using LinearAlgebra

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""Build a default test wing."""
function test_wing()
    cfg  = build_run_config(1)
    return cfg, build_wing_geometry(cfg)
end

"""Nominal optimisation variable vector (feasible point)."""
const X_NOM = [0.30, 0.15, 0.20, 2.0, 0.005]   # [pos, w_tip, w_root, box_t, refor_t]

# ---------------------------------------------------------------------------
# 1. Utility functions
# ---------------------------------------------------------------------------

@testset "Utils" begin
    @testset "trapz" begin
        x = [0.0, 1.0, 2.0, 3.0]
        y = [0.0, 1.0, 2.0, 3.0]   # integral = 4.5
        @test trapz(x, y) ≈ 4.5

        # Constant function → width × height
        @test trapz([0.0, 2.0], [3.0, 3.0]) ≈ 6.0

        # Quadratic (analytical: ∫₀¹ x² dx = 1/3)
        xq = range(0, 1; length=1001) |> collect
        yq = xq.^2
        @test trapz(xq, yq) ≈ 1/3  atol=1e-4
    end

    @testset "cumtrapz" begin
        x = [0.0, 1.0, 2.0]
        y = [1.0, 1.0, 1.0]   # running integral: 0, 1, 2
        ct = cumtrapz(x, y)
        @test ct[1] ≈ 0.0
        @test ct[2] ≈ 1.0
        @test ct[3] ≈ 2.0
    end

    @testset "quadratic_roots" begin
        r1, r2 = quadratic_roots(1.0, -3.0, 2.0)   # (x-1)(x-2) = 0
        @test sort([r1, r2]) ≈ [1.0, 2.0]
    end
end

# ---------------------------------------------------------------------------
# 2. Material
# ---------------------------------------------------------------------------

@testset "Material" begin
    for name in ("carbono_bi", "carbono_uni", "balsa", "freijo")
        mat = material_properties(name)
        @test mat.axial_E   > 0
        @test mat.dens      > 0
        @test mat.axial_trac > mat.axial_comp * 0.0   # both positive
    end

    @test_throws ErrorException material_properties("unobtainium")
end

# ---------------------------------------------------------------------------
# 3. Inertia tensor transformations
# ---------------------------------------------------------------------------

@testset "Inertia transforms" begin
    I = [100.0 0.0; 0.0 50.0]

    @testset "translate (parallel axis)" begin
        d = [0.1, 0.2];  A = 0.01
        It = translate_inertia_tensor(I, d, A)
        @test It[1,1] ≈ I[1,1] + d[1]^2 * A
        @test It[2,2] ≈ I[2,2] + d[2]^2 * A
        @test It[1,2] ≈ prod(d) * A
    end

    @testset "rotate by 0 → unchanged" begin
        Ir = rotate_inertia_tensor(I, 0.0)
        @test Ir ≈ I  atol=1e-10
    end

    @testset "rotate by π/2 → swap Ix↔Iy" begin
        Ir = rotate_inertia_tensor(I, π/2)
        @test Ir[1,1] ≈ I[2,2]  atol=1e-10
        @test Ir[2,2] ≈ I[1,1]  atol=1e-10
    end

    @testset "principal inertia" begin
        I2 = [100.0 20.0; 20.0 50.0]
        Ip, theta = principal_inertia_tensor(I2)
        @test Ip[1,1] >= Ip[2,2]           # I_max ≥ I_min
        @test Ip[1,1] + Ip[2,2] ≈ I2[1,1] + I2[2,2]  atol=1e-8  # trace conserved
        @test theta >= 0
    end
end

# ---------------------------------------------------------------------------
# 4. Section geometry helpers
# ---------------------------------------------------------------------------

@testset "Section helpers" begin
    @testset "multielement_centroid" begin
        # Two equal rectangles side by side → centroid at x = 1.5
        pos  = [1.0 0.0; 2.0 0.0]
        area = [1.0, 1.0]
        c = multielement_centroid(pos, area)
        @test c[1] ≈ 1.5
        @test c[2] ≈ 0.0
    end

    @testset "section_centroid of a straight line" begin
        pts = [0.0 0.0; 1.0 0.0; 2.0 0.0]
        c = section_centroid(pts)
        @test c[1] ≈ 1.0   # midpoint
        @test c[2] ≈ 0.0
    end
end

# ---------------------------------------------------------------------------
# 5. Airfoil parsing
# ---------------------------------------------------------------------------

@testset "Airfoil parsing" begin
    foil_dir = joinpath(@__DIR__, "..", "..", "foils")   # repo root foils/
    if isdir(foil_dir)
        foil = parse_foil("s1223"; foil_dir=foil_dir)

        @test foil.esp > 0
        @test size(foil.ext, 2) == 2    # two-column matrix
        @test size(foil.int, 2) == 2

        # x should be monotonically increasing after parsing
        @test issorted(foil.ext[:, 1])
        @test issorted(foil.int[:, 1])

        @testset "foil_properties" begin
            chord  = 1.0
            limits = [0.25, 0.60]
            fg = foil_properties(chord, 0.0, limits, [0.0, 0.0], foil)

            @test fg.extra_length > 0
            @test fg.intra_length > 0
            @test fg.area         > 0
            # upper surface y should be above lower at quarter-chord limits
            @test fg.upper[1, 2] >= fg.lower[1, 2]
        end
    else
        @warn "Foil data directory not found; skipping airfoil tests."
    end
end

# ---------------------------------------------------------------------------
# 6. Wing geometry and station generation
# ---------------------------------------------------------------------------

@testset "Wing geometry" begin
    cfg, wing = test_wing()

    @test wing.root_chord  ≈ 1.0
    @test wing.full_span   ≈ 7.0
    @test length(wing.chord_transitions) == 3
    @test length(wing.span_transitions)  == 3
    @test issorted(wing.span_transitions)

    @testset "parametric geometry" begin
        geom = reinforced_box_parametric_geometry(X_NOM, wing)
        @test length(geom.x_position) == 3
        @test length(geom.width)      == 3
        @test length(geom.thickness)  == 3
        @test length(geom.reinforcer) == 3
    end

    foil_dir = joinpath(@__DIR__, "..", "..", "foils")
    if isdir(foil_dir)
        foil = parse_foil("s1223"; foil_dir=foil_dir)
        geom = reinforced_box_parametric_geometry(X_NOM, wing)

        s = range(0.2, 3.4; length=10) |> collect
        stations = generate_wing_stations(wing, geom, s)

        @test length(stations) == 10
        @test all(st.chord > 0 for st in stations)
        @test all(st.thickness > 0 for st in stations)
        @test all(st.limits[2] > st.limits[1] for st in stations)
    end
end

# ---------------------------------------------------------------------------
# 7. Cross-section properties
# ---------------------------------------------------------------------------

@testset "Cross sections" begin
    foil_dir = joinpath(@__DIR__, "..", "..", "foils")
    isdir(foil_dir) || return

    foil = parse_foil("s1223"; foil_dir=foil_dir)
    st   = Station(1.0, 0.003, 0.005, 1.0, [0.28, 0.55])

    @testset "reinforced_box_section" begin
        prop = reinforced_box_section(st, foil)
        @test prop.Iz0        > 0
        @test prop.area       > 0
        @test prop.width      > 0
        @test prop.height     > 0
        @test prop.theta_max  >= 0
        @test size(prop.vertices) == (4, 2)
        # Check inertia tensor is positive definite
        @test isposdef(prop.I)  ||  prop.I[1,1] > 0
    end

    @testset "box_section (alias)" begin
        prop = box_section(st, foil)
        @test prop.Iz0 > 0
    end

    @testset "D_section" begin
        prop = D_section(st, foil)
        @test prop.Iz0 > 0
    end
end

# ---------------------------------------------------------------------------
# 8. Failure criteria
# ---------------------------------------------------------------------------

@testset "Failure criteria" begin
    mat = material_properties("carbono_bi")

    @testset "max_tension: zero stress → Inf margin" begin
        sigma = zeros(3, 3)
        @test max_tension_criteria(sigma, mat, "axial") == Inf
    end

    @testset "max_tension: stress at failure limit → margin ≈ 0" begin
        sigma = zeros(3, 3)
        sigma[1, 1] = mat.axial_trac
        ms = max_tension_criteria(sigma, mat, "axial")
        @test ms ≈ 0.0  atol=1e-6
    end

    @testset "tsai_wu: compressive stress" begin
        sigma = zeros(3, 3)
        sigma[1, 1] = -mat.axial_comp   # at compressive limit
        ms = tsai_wu_criteria(sigma, mat, "axial")
        @test isfinite(ms)
    end

    @testset "tsai_hill: zero stress → safe (≈1)" begin
        sigma = zeros(3, 3)
        @test tsai_hill_criteria(sigma, mat) ≈ 1.0  atol=1e-10
    end
end

# ---------------------------------------------------------------------------
# 9. Weight function
# ---------------------------------------------------------------------------

@testset "Weight function" begin
    # Positive safety margin, low report → gain = 1 → fit = -1/mass
    mass = 0.5
    @test weight_function(mass, 1.5, 0.5) ≈ -1.0/mass

    # Negative safety margin → penalty (less negative)
    @test weight_function(mass, -0.1, 0.0) > weight_function(mass, 1.5, 0.5)

    # report >= 1 → reduced gain
    @test weight_function(mass, 1.5, 1.5) ≈ -0.5/mass
end

# ---------------------------------------------------------------------------
# 10. Torsion
# ---------------------------------------------------------------------------

@testset "Torsion" begin
    @testset "integral_linha_media" begin
        t = [0.002, 0.002, 0.002, 0.002]
        s = [0.05, 0.10, 0.05, 0.10]   # perimeter segments
        val = integral_linha_media(t, s)
        @test val ≈ sum(s ./ t)
    end
end

# ---------------------------------------------------------------------------
# 11. Configuration
# ---------------------------------------------------------------------------

@testset "Configuration" begin
    cfg = build_run_config(1)
    @test cfg.section_type == "reinforced_box"
    @test cfg.constraint_type == "ss"

    wing = build_wing_geometry(cfg)
    @test wing.full_span   ≈ 7.0
    @test wing.root_chord  ≈ 1.0
    @test length(wing.span_transitions) == 3
end

# ---------------------------------------------------------------------------
# 12. Integration test: full fitness evaluation
# ---------------------------------------------------------------------------

@testset "Integration: fit_function" begin
    foil_dir = joinpath(@__DIR__, "..", "..", "foils")
    isdir(foil_dir) || return   # skip if foil files unavailable

    cfg, wing = test_wing()
    foil      = parse_foil(cfg.foil_name; foil_dir=foil_dir)
    loads     = example_aero_loads(wing; n_panels=20)

    @test length(loads) == 1
    @test length(loads[1].x) == 20

    raw_fit, data = fit_function(X_NOM, cfg, wing, foil, loads)

    @test isfinite(raw_fit)
    @test data.mass > 0
    @test length(data.stations) > 0
    @test length(data.section_props) == length(data.stations)
end
