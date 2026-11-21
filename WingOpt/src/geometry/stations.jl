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
    length(geometry.x_position) == length(wing.span_transitions) || throw(
        DimensionMismatch("spar parameters must be supplied at every wing transition"),
    )
    length(s) > 0 || return Station[]
    all(wing.span_transitions[1] .<= s .<= wing.span_transitions[end]) ||
        throw(ArgumentError("station spans must lie within the modeled semi-span"))
    function interp(values, span)
        i = clamp(
            searchsortedlast(wing.span_transitions, span),
            1,
            length(wing.span_transitions) - 1,
        )
        x1, x2 = wing.span_transitions[i], wing.span_transitions[i+1]
        t = (span - x1) / (x2 - x1)
        return muladd(t, values[i+1] - values[i], values[i])
    end
    return [
        begin
            span = Float64(si)
            chord = interp(wing.chord_transitions, span)
            front = (interp(geometry.x_position, span) + chord / 4) / chord
            width = interp(geometry.width, span)
            Station(
                chord,
                interp(geometry.thickness, span),
                interp(geometry.reinforcer, span),
                span,
                [front, front + width / chord],
            )
        end for si in s
    ]
end
