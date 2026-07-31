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
function generate_wing_stations(wing::WingGeometry,
                                 geometry::StructuredGeometry,
                                 s::AbstractVector{<:Real})
    sp = wing.span_transitions       # control-point spans
    ch = wing.chord_transitions      # control-point chords

    # Local chord at each station
    chord_itp = LinearInterpolation(sp, ch)
    chord = chord_itp.(s)

    # Front-web position → fractional chord limit
    #   x_position is in metres, offset from c/4 (same convention as basic_box_section)
    xpos_itp = LinearInterpolation(sp, geometry.x_position)
    a_limit_m = xpos_itp.(s)
    a_limit   = (a_limit_m .+ chord ./ 4) ./ chord   # fractional [0,1]

    # Width → back fractional limit
    width_itp = LinearInterpolation(sp, geometry.width)
    w_m    = width_itp.(s)
    b_limit = a_limit .+ w_m ./ chord

    # Thicknesses
    thick_itp = LinearInterpolation(sp, geometry.thickness)
    t_box     = thick_itp.(s)

    refor_itp = LinearInterpolation(sp, geometry.reinforcer)
    t_refor   = refor_itp.(s)

    return [Station(chord[i], t_box[i], t_refor[i], Float64(s[i]),
                    [a_limit[i], b_limit[i]])
            for i in eachindex(s)]
end
