"""
Wing station generation.
Translates `geometry/generate_wing_stations.m`.
Uses Interpolations.LinearInterpolation instead of Matlab's interp1.
"""

"""
    generate_wing_stations(wing, geometry, s) -> Vector{Station}

Interpolate cross-section parameters at spanwise positions `s` and
return one `Station` per position.

`wing`     – WingGeometry (provides span_transitions and chord_transitions)
`geometry` – StructuredGeometry (spar parametrisation, one value per span break)
`s`        – vector of desired spanwise positions (m, right semi-span)
"""
function generate_wing_stations(
    wing::WingGeometry,
    geometry::StructuredGeometry,
    s::AbstractVector{<:Real},
)
    transitions = wing.span_transitions
    n = length(transitions)
    n >= 2 || throw(ArgumentError("at least two span transitions are required"))
    length(wing.chord_transitions) == n ||
        throw(DimensionMismatch("chord and span transitions must have equal lengths"))
    all(
        length(values) == n for values in
        (geometry.x_position, geometry.width, geometry.thickness, geometry.reinforcer)
    ) || throw(
        DimensionMismatch("all structural geometry arrays must match span transitions"),
    )
    all(diff(transitions) .> 0) ||
        throw(ArgumentError("span transitions must be strictly increasing"))
    isempty(s) && return Station[]
    all(isfinite, s) || throw(ArgumentError("station positions must be finite"))
    all((transitions[1] .<= s) .& (s .<= transitions[end])) ||
        throw(DomainError(s, "station positions must lie within the wing semi-span"))

    chord = _linear_sample(transitions, wing.chord_transitions, s)
    front_x = _linear_sample(transitions, geometry.x_position, s)
    width = _linear_sample(transitions, geometry.width, s)
    thickness = _linear_sample(transitions, geometry.thickness, s)
    reinforcer = _linear_sample(transitions, geometry.reinforcer, s)

    stations = Vector{Station}(undef, length(s))
    for i in eachindex(s)
        chord[i] > 0 || throw(DomainError(chord[i], "local chord must be positive"))
        front = (front_x[i] + chord[i] / 4) / chord[i]
        back = front + width[i] / chord[i]
        0 <= front < back <= 1 || throw(
            DomainError((front, back), "spar limits must satisfy 0 ≤ front < back ≤ 1"),
        )
        thickness[i] > 0 ||
            throw(DomainError(thickness[i], "wall thickness must be positive"))
        reinforcer[i] >= 0 ||
            throw(DomainError(reinforcer[i], "reinforcer thickness must be nonnegative"))
        stations[i] =
            Station(chord[i], thickness[i], reinforcer[i], Float64(s[i]), [front, back])
    end
    return stations
end

function _linear_sample(x::AbstractVector, y::AbstractVector, points::AbstractVector)
    length(x) == length(y) ||
        throw(DimensionMismatch("interpolation knots and values must match"))
    return [
        begin
            i = min(searchsortedlast(x, point), length(x) - 1)
            fraction = (point - x[i]) / (x[i+1] - x[i])
            Float64(y[i] + fraction * (y[i+1] - y[i]))
        end for point in points
    ]
end
