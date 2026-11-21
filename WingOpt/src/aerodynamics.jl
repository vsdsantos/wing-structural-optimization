module Aerodynamics

using ..WingOpt: AeroLoad

export elliptic_load

"""
    elliptic_load(semispan, total_lift, n_panels) -> AeroLoad

Create a symmetric elliptic lift distribution over a semi-span. `total_lift`
is the combined lift on both wing halves.
"""
function elliptic_load(semispan::Real, total_lift::Real, n_panels::Integer)
    isfinite(semispan) && semispan > 0 ||
        throw(ArgumentError("semispan must be finite and positive"))
    isfinite(total_lift) && total_lift >= 0 ||
        throw(ArgumentError("total_lift must be finite and nonnegative"))
    n_panels > 0 || throw(ArgumentError("n_panels must be positive"))

    span = Float64(semispan)
    panel_width = span / n_panels
    x = [(i - 0.5) * panel_width for i in 1:n_panels]
    shape = sqrt.(1 .- (x ./ span) .^ 2)
    normalization = 2 * sum(shape) * panel_width
    lift_per_span = Float64(total_lift) .* shape ./ normalization

    return AeroLoad(x, lift_per_span, zeros(n_panels))
end

end
