"""
Airfoil parsing and geometric property evaluation.
Replaces `geometry/foil/parse_foil.m` and `foil_properties.m`.

Key differences from the Matlab version:
- Uses Interpolations.jl LinearInterpolation instead of the Curve Fitting Toolbox.
- Returns in-memory FoilData instead of writing a .mat cache.
- `foil_properties` accepts a FoilData struct instead of a file name.
"""

"""
    parse_foil(foil_name; foil_dir="foils") -> FoilData

Parse a Selig-format `.dat` airfoil file and build linear interpolants for
the upper (extrados) and lower (intrados) surfaces.

The file is expected to live at `foil_dir/foil_name.dat`.
The first non-numeric line (title) is skipped automatically.
"""
function parse_foil(foil_name::String; foil_dir::String="foils")
    error("Not implemented")
end

"""
    foil_properties(chord, scale_thick, limits, offset, foil) -> FoilGeometry

Evaluate airfoil geometry at a spanwise station.

Arguments:
- `chord`       – local chord (m)
- `scale_thick` – thickness scaling offset (set to 0 for normal use)
- `limits`      – [a, b] front/back limits in fractional chord [0, 1]
- `offset`      – [upper_offset, lower_offset] y offsets (use [0,0] normally)
- `foil`        – FoilData from `parse_foil`

Returns a FoilGeometry with upper/lower surface arrays in the quarter-chord
coordinate system (x = 0 at c/4, positive toward TE; y = 0 at camber line).
"""
function foil_properties(chord::Float64, scale_thick::Float64,
                          limits::AbstractVector, offset::AbstractVector,
                          foil::FoilData)
    error("Not implemented")
end
