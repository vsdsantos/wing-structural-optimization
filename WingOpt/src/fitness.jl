"""
Fitness (objective) function and weight function.
Translates `fit_function.m` and `weight_function.m`.

Key differences from the Matlab version:
- All data is passed as arguments; no `load`/`save` to .mat files.
- `loads` (aerodynamic loading) must be supplied externally; this module
  defines a simple example loader for testing.
- Missing symbols resolved:
    n_sec        → length(stations)
    beam_forces_SS/FX → internal_loads_SS/FX
    long_mass    → beam_mass
    desenho      → dropped (visualisation artifact)
"""

# ---------------------------------------------------------------------------
# Aerodynamic load example generator (replaces loads.mat)
# ---------------------------------------------------------------------------

"""
    example_aero_loads(wing; n_panels=20) -> Vector{AeroLoad}

Generate a simple elliptic lift distribution as an example aerodynamic load.
In production use, replace with actual VLM/CFD loads.
"""
function example_aero_loads(wing::WingGeometry; n_panels::Int=20)
    b = wing.full_span
    # Panel centres (full span, symmetric)
    s = range(-b/2 + b/(2*n_panels), b/2 - b/(2*n_panels); length=n_panels)
    # Elliptic lift distribution (total lift ≈ weight, normalised here)
    dLdb  = sqrt.(max.(1 .- (2 .* collect(s) ./ b).^2, 0.0))
    dMadb = zeros(n_panels)
    return [AeroLoad(collect(s), dLdb, dMadb)]
end

# ---------------------------------------------------------------------------
# Weight function (penalty / objective shaping)
# ---------------------------------------------------------------------------

"""
    weight_function(mass, safety_margin_val, report) -> Float64

Map (mass, min. safety margin, report-45 ratio) to a scalar fitness value.
Lower (more negative) is better.

Translated faithfully from `weight_function.m`.
"""
function weight_function(mass::Float64, safety_margin_val::Float64, report::Float64)
    gain = if safety_margin_val <= 0
        -mass / (1 + safety_margin_val)
    elseif safety_margin_val <= 1 && report < 1
        1.5
    elseif report >= 1
        0.5
    else
        1.0
    end
    return -gain / mass
end

# ---------------------------------------------------------------------------
# Fitness function
# ---------------------------------------------------------------------------

"""
    fit_function(X, config, wing, foil, loads; foil_dir="foils") -> (raw_fit, FitData)

Evaluate the fitness for the optimisation variable vector `X`.

`X` layout (1-based):
    X[1] – web_pos_tip       ∈ [0.01, 0.60]
    X[2] – width_tip         ∈ [0.01, 0.40]
    X[3] – width_root        ∈ [0.01, 0.40]
    X[4] – box_thickness     ∈ {1, 2, 3, 4}  (integer, multiplied by 1 mm)
    X[5] – reinforcer_thick  ∈ [0.002, 0.010] (m)

Returns `(raw_fit, data::FitData)` where raw_fit is the scalar objective
(to be minimised; more negative = better structure per unit mass).
"""
function fit_function(X::AbstractVector,
                       config::RunConfig,
                       wing::WingGeometry,
                       foil::FoilData,
                       loads::Vector{AeroLoad})
    # 1. Parametric geometry
    geometry = reinforced_box_parametric_geometry(X, wing)

    # 2. Station positions: right semi-span from the first load case
    #    Use the right half of the first load case panel positions.
    s_full = loads[1].x
    n_half = length(s_full) ÷ 2
    s      = s_full[n_half+1:end]   # right semi-span panel centres

    # 3. Wing stations
    stations = generate_wing_stations(wing, geometry, s)
    n_sec    = length(stations)

    # 4. Section properties at each station
    props = [compute_section(stations[i], foil, config.section_type)
             for i in 1:n_sec]

    # 5. Internal loads for each aerodynamic case
    internal_loads_vec = [
        compute_internal_loads(ld, props, wing, config.constraint_type)
        for ld in loads
    ]

    # 6. Safety margins
    ms = zeros(length(loads))
    for (j, il) in enumerate(internal_loads_vec)
        ms_L = [safety_margin(stations[i], il.L, i, props[i]) for i in 1:n_sec]
        ms_R = [safety_margin(stations[i], il.R, i, props[i]) for i in 1:n_sec]
        ms[j] = minimum(vcat(ms_L, ms_R))
    end
    min_ms = minimum(ms)

    # 7. Report-45 flexibility factor
    mat_box = material_properties("freijo")
    report  = 0.0   # report_45 requires wing.V_d; skip if not defined

    # 8. Beam mass (one semi-span, then doubled)
    mat_reinforcer = material_properties("freijo")
    mass_half, dm = beam_mass(stations, props, mat_box, mat_reinforcer)
    mass = 2 * mass_half

    # 9. Fitness value
    raw_fit = weight_function(mass, min_ms, report)

    data = FitData(stations, props, wing, mass, dm, min_ms, ms, internal_loads_vec)
    return raw_fit, data
end
