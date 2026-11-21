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
The lift and pitching-moment distributions are normalized to `total_lift`
(N) and `total_pitching_moment` (N·m), respectively. Panel positions are
full-span centers; an even `n_panels` is required for a symmetric wing.
In production use, replace with actual VLM/CFD loads.
"""
function example_aero_loads(
    wing::WingGeometry;
    n_panels::Int = 20,
    total_lift::Real = 1.0,
    total_pitching_moment::Real = 0.0,
)
    n_panels > 0 && iseven(n_panels) ||
        throw(ArgumentError("n_panels must be a positive even number"))
    isfinite(wing.full_span) && wing.full_span > 0 ||
        throw(ArgumentError("wing full span must be positive and finite"))
    isfinite(total_lift) && total_lift >= 0 ||
        throw(ArgumentError("total_lift must be nonnegative and finite"))
    isfinite(total_pitching_moment) ||
        throw(ArgumentError("total_pitching_moment must be finite"))

    panel_width = wing.full_span / n_panels
    x = collect(
        range(-wing.full_span / 2 + panel_width / 2; step = panel_width, length = n_panels),
    )
    shape = sqrt.(max.(0.0, 1 .- (2 .* x ./ wing.full_span) .^ 2))
    normalization =
        trapz(vcat(-wing.full_span / 2, x, wing.full_span / 2), vcat(0.0, shape, 0.0))
    dLdb = shape .* (total_lift / normalization)
    dMadb = shape .* (total_pitching_moment / normalization)
    return [AeroLoad(x, dLdb, dMadb)]
end

"""
    weight_function(mass, safety_margin_val, report) -> Float64

Map (mass, min. safety margin, report-45 ratio) to a scalar fitness value.
Lower (more negative) is better.

Translated faithfully from `weight_function.m`.
"""
function weight_function(mass::Float64, safety_margin_val::Float64, report::Float64)
    error("Not implemented")
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
    error("Not implemented")
end
