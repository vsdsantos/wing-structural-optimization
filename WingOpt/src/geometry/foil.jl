"""
Airfoil parsing and geometric property evaluation.
Replaces `geometry/foil/parse_foil.m` and `foil_properties.m`.

Key differences from the Matlab version:
- Uses Interpolations.jl gridded linear interpolation instead of the Curve Fitting Toolbox.
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
function parse_foil(foil_name::String; foil_dir::String = "foils")
    path = joinpath(foil_dir, endswith(foil_name, ".dat") ? foil_name : "$foil_name.dat")
    isfile(path) || throw(ArgumentError("airfoil file not found: $path"))
    points = Vector{Vector{Float64}}()
    for line in eachline(path)
        values = tryparse.(Float64, split(strip(line)))
        length(values) >= 2 && all(value -> !isnothing(value), values[1:2]) &&
            push!(points, [values[1]::Float64, values[2]::Float64])
    end
    length(points) >= 4 ||
        throw(ArgumentError("airfoil file contains too few coordinate points"))
    data = reduce(vcat, permutedims.(points))
    leading_edge = argmin(data[:, 1])
    leading_edge > 1 && leading_edge < size(data, 1) ||
        throw(ArgumentError("airfoil coordinates must contain upper and lower surfaces"))
    upper = _unique_foil_points(data[1:leading_edge, :])
    lower = _unique_foil_points(data[leading_edge:end, :])
    size(upper, 1) >= 2 && size(lower, 1) >= 2 ||
        throw(
            ArgumentError(
                "each airfoil surface must contain at least two unique x coordinates",
            ),
        )
    spacing = min(minimum(diff(upper[:, 1])), minimum(diff(lower[:, 1]))) / 10
    return FoilData(
        interpolate((upper[:, 1],), upper[:, 2], Gridded(Linear())),
        interpolate((lower[:, 1],), lower[:, 2], Gridded(Linear())),
        upper,
        lower,
        spacing,
    )
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
function foil_properties(
    chord::Float64,
    scale_thick::Float64,
    limits::AbstractVector,
    offset::AbstractVector,
    foil::FoilData,
)
    chord > 0 || throw(DomainError(chord, "chord must be positive"))
    length(limits) == 2 && length(offset) == 2 ||
        throw(DimensionMismatch("limits and offset must each have two values"))
    0 <= limits[1] < limits[2] <= 1 ||
        throw(DomainError(limits, "limits must be increasing fractions within [0, 1]"))
    foil.esp > 0 ||
        throw(DomainError(foil.esp, "airfoil evaluation spacing must be positive"))

    reference_x = collect(0.0:foil.esp:1.0)
    isempty(reference_x) && (reference_x = [0.0, 1.0])
    shift_index = argmin(abs.(reference_x .- 0.25))
    shift = (
        foil.ex_lin(reference_x[shift_index]) + foil.in_lin(reference_x[shift_index])
    ) / 2
    x_fraction = collect(limits[1]:foil.esp:limits[2])
    isempty(x_fraction) && (x_fraction = [limits[1]])
    if !isapprox(last(x_fraction), limits[2]; atol = eps(Float64) * 8, rtol = 0)
        push!(x_fraction, limits[2])
    end
    upper = hcat(
        x_fraction .* chord .- chord / 4,
        [foil.ex_lin(x) * chord - shift * chord for x in x_fraction],
    )
    lower = hcat(
        x_fraction .* chord .- chord / 4,
        [foil.in_lin(x) * chord - shift * chord for x in x_fraction],
    )
    if !iszero(scale_thick)
        upper = _offset_foil_surface(upper, scale_thick)
        lower = _offset_foil_surface(lower, -scale_thick)
    end
    upper[:, 2] .-= offset[1]
    lower[:, 2] .+= offset[2]

    upper_length = sum(hypot.(diff(upper[:, 1]), diff(upper[:, 2])))
    lower_length = sum(hypot.(diff(lower[:, 1]), diff(lower[:, 2])))
    front_length = upper[1, 2] - lower[1, 2]
    back_length = upper[end, 2] - lower[end, 2]
    area = trapz(lower[:, 1], upper[:, 2] .- lower[:, 2])
    x = vcat(
        upper[:, 1], lower[:, 1], range(lower[1, 1], upper[1, 1]; length = 20),
        range(lower[end, 1], upper[end, 1]; length = 20),
    )
    y = vcat(
        upper[:, 2], lower[:, 2], range(lower[1, 2], upper[1, 2]; length = 20),
        range(lower[end, 2], upper[end, 2]; length = 20),
    )
    return FoilGeometry(
        upper, lower, upper_length, lower_length, front_length, back_length,
        upper_length + lower_length + front_length + back_length, area, x, y,
    )
end

function _unique_foil_points(points::AbstractMatrix)
    ordered = points[sortperm(points[:, 1]), :]
    keep = [true; diff(ordered[:, 1]) .> eps(Float64) * 8]
    return ordered[keep, :]
end

function _offset_foil_surface(points::Matrix{Float64}, distance::Float64)
    size(points, 1) >= 2 || return points
    tangent = diff(points; dims = 1)
    lengths = hypot.(tangent[:, 1], tangent[:, 2])
    all(>(0), lengths) ||
        throw(ArgumentError("airfoil surface contains duplicate adjacent points"))
    midpoint = points[1:(end - 1), :] .+ tangent ./ 2
    normal = hcat(-tangent[:, 2], tangent[:, 1]) ./ lengths
    return midpoint .+ distance .* normal
end
