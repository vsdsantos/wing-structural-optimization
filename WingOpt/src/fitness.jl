"""
Fitness (objective) function and weight function.
Translates `fit_function.m` and `weight_function.m`.

Key differences from the Matlab version:
- All data is passed as arguments; no `load`/`save` to .mat files.
- `loads` (aerodynamic loading) must be supplied externally; this module
  defines a simple example loader for testing.
"""

"""
    example_aero_loads(wing; n_panels=20) -> Vector{AeroLoad}

Generate a simple elliptic lift distribution as an example aerodynamic load.
In production use, replace with actual VLM/CFD loads.
"""
function example_aero_loads(wing::WingGeometry; n_panels::Int = 20)
    n_panels >= 2 || throw(ArgumentError("n_panels must be at least two"))
    wing.full_span > 0 || throw(ArgumentError("wing span must be positive"))
    step = wing.full_span / n_panels
    x = collect(
        range(
            -wing.full_span / 2 + step / 2,
            wing.full_span / 2 - step / 2;
            length = n_panels,
        ),
    )
    shape = sqrt.(max.(0.0, 1 .- (2 .* x ./ wing.full_span) .^ 2))
    q_peak = 4 * 100.0 / (pi * wing.full_span)
    lift = q_peak .* shape
    return [AeroLoad(x, lift, zeros(n_panels))]
end

"""
    weight_function(mass, safety_margin_val, report) -> Float64

Map (mass, min. safety margin, report-45 ratio) to a scalar fitness value.
Lower (more negative) is better.

Translated faithfully from `weight_function.m`.
"""
function weight_function(mass::Float64, safety_margin_val::Float64, report::Float64)
    mass > 0 && isfinite(mass) || throw(ArgumentError("mass must be positive and finite"))
    if safety_margin_val <= 0
        gain = -mass / (1 + safety_margin_val)
    elseif report >= 1
        gain = 0.5
    elseif safety_margin_val <= 1
        gain = 1.5
    else
        gain = 1.0
    end
    return -gain / mass
end

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
function fit_function(
    X::AbstractVector,
    config::RunConfig,
    wing::WingGeometry,
    foil::FoilData,
    loads::Vector{AeroLoad},
)
    length(X) == 5 || throw(DimensionMismatch("design vector X must contain five values"))
    isempty(loads) && throw(ArgumentError("at least one aerodynamic load case is required"))
    geometry = reinforced_box_parametric_geometry(X, wing)
    span = sort(unique(reduce(vcat, (load.x[load.x .>= 0] for load in loads))))
    isempty(span) && throw(ArgumentError("loads must include right-wing spanwise stations"))
    span[end] <= wing.full_span / 2 ||
        throw(ArgumentError("load stations exceed the semi-span"))
    stations = generate_wing_stations(wing, geometry, span)
    section_props =
        [compute_section(station, foil, config.section_type) for station in stations]
    box_material = material_properties("balsa")
    reinforcer_material = material_properties("freijo")
    mass_half, dm_half =
        beam_mass(stations, section_props, box_material, reinforcer_material)
    mass = 2mass_half
    internal_loads = [
        compute_internal_loads(load, section_props, wing, config.constraint_type) for
        load in loads
    ]
    ms = [
        minimum(
            vcat(
                [
                    safety_margin(stations[i], case.R, i, section_props[i]) for
                    i in eachindex(stations)
                ],
                [
                    safety_margin(stations[i], case.L, i, section_props[i]) for
                    i in eachindex(stations)
                ],
            ),
        ) for case in internal_loads
    ]
    min_ms = minimum(ms)
    report = report_45(stations, section_props, wing, box_material)
    data = FitData(
        stations,
        section_props,
        wing,
        mass,
        2 .* dm_half,
        min_ms,
        ms,
        internal_loads,
    )
    return weight_function(mass, min_ms, report), data
end
