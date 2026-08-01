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
    error("Not implemented")
end
