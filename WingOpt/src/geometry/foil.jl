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
function parse_foil(foil_name::String; foil_dir::String = "foils")
    isempty(strip(foil_name)) && throw(ArgumentError("airfoil name must not be empty"))
    path = joinpath(
        foil_dir,
        endswith(lowercase(foil_name), ".dat") ? foil_name : "$foil_name.dat",
    )
    if !isfile(path)
        packaged_foils = normpath(joinpath(@__DIR__, "../../../foils"))
        candidate = joinpath(
            packaged_foils,
            endswith(lowercase(foil_name), ".dat") ? foil_name : "$foil_name.dat",
        )
        isfile(candidate) && (path = candidate)
    end
    isfile(path) || throw(ArgumentError("airfoil file not found: $path"))

    coordinates = NTuple{2,Float64}[]
    for line in eachline(path)
        fields = split(strip(line))
        length(fields) >= 2 || continue
        x, y = tryparse(Float64, fields[1]), tryparse(Float64, fields[2])
        if x !== nothing && y !== nothing
            push!(coordinates, (x, y))
        end
    end
    length(coordinates) >= 5 ||
        throw(ArgumentError("airfoil file must contain at least five coordinate pairs"))
    data = collect.(coordinates)
    all(p -> all(isfinite, p), data) ||
        throw(ArgumentError("airfoil coordinates must be finite"))
    all(p -> 0 <= p[1] <= 1, data) ||
        throw(DomainError(data, "airfoil x-coordinates must lie in [0, 1]"))

    leading_edge = argmin(first.(coordinates))
    leading_edge > 1 && leading_edge < length(coordinates) || throw(
        ArgumentError(
            "airfoil coordinates must trace upper and lower surfaces from a trailing edge",
        ),
    )
    upper = _sorted_surface(data[1:leading_edge])
    lower = _sorted_surface(data[leading_edge:end])
    spacing = min(minimum(diff(upper[:, 1])), minimum(diff(lower[:, 1])))
    spacing > 0 || throw(
        ArgumentError("airfoil coordinates must contain distinct x-values on each surface"),
    )
    return FoilData(
        linear_interpolation(upper[:, 1], upper[:, 2]; extrapolation_bc = Line()),
        linear_interpolation(lower[:, 1], lower[:, 2]; extrapolation_bc = Line()),
        upper,
        lower,
        spacing / 10,
    )
end

"""
    foil_properties(chord, scale_thick, limits, offset, foil) -> FoilGeometry

Evaluate airfoil geometry at a spanwise station.

Arguments:
- `chord`       – local chord (m)
- `scale_thick` – signed normal offset (m; positive moves both surfaces outward)
- `limits`      – [a, b] front/back limits in fractional chord [0, 1]
- `offset`      – [upper, lower] y offsets (m; upper is subtracted, lower added)
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
    isfinite(chord) && chord > 0 ||
        throw(DomainError(chord, "chord must be positive and finite"))
    isfinite(scale_thick) || throw(ArgumentError("thickness scale must be finite"))
    length(limits) == 2 || throw(DimensionMismatch("limits must be [front, back]"))
    length(offset) == 2 || throw(DimensionMismatch("offset must be [upper, lower]"))
    all(isfinite, limits) && all(isfinite, offset) ||
        throw(ArgumentError("limits and offsets must be finite"))
    0 <= limits[1] < limits[2] <= 1 ||
        throw(DomainError(limits, "limits must satisfy 0 ≤ front < back ≤ 1"))

    count = max(2, ceil(Int, (limits[2] - limits[1]) / foil.esp) + 1)
    x_chord = collect(range(Float64(limits[1]), Float64(limits[2]); length = count))
    common_min = max(first(foil.ext[:, 1]), first(foil.int[:, 1]))
    common_max = min(last(foil.ext[:, 1]), last(foil.int[:, 1]))
    x_sample = clamp.(x_chord, common_min, common_max)
    shift = (foil.ex_lin(0.25) + foil.in_lin(0.25)) / 2
    upper = hcat(
        chord .* (x_chord .- 0.25),
        chord .* ([foil.ex_lin(x) for x in x_sample] .- shift),
    )
    lower = hcat(
        chord .* (x_chord .- 0.25),
        chord .* ([foil.in_lin(x) for x in x_sample] .- shift),
    )

    if !iszero(scale_thick)
        _offset_surface!(upper, scale_thick, 1)
        _offset_surface!(lower, scale_thick, -1)
    end
    upper[:, 2] .-= offset[1]
    lower[:, 2] .+= offset[2]
    gap = upper[:, 2] - lower[:, 2]
    gap_tolerance = 1e-12 * chord
    (all(gap .>= -gap_tolerance) && any(gap .> gap_tolerance)) ||
        throw(DomainError(offset, "surface offsets/thickness must not invert the airfoil"))

    extra_length = sum(hypot.(diff(upper[:, 1]), diff(upper[:, 2])))
    intra_length = sum(hypot.(diff(lower[:, 1]), diff(lower[:, 2])))
    front_length = upper[1, 2] - lower[1, 2]
    back_length = upper[end, 2] - lower[end, 2]
    area = trapz(upper[:, 1], upper[:, 2] - lower[:, 2])
    perimeter = extra_length + intra_length + front_length + back_length
    x = vcat(
        upper[:, 1],
        lower[:, 1],
        range(lower[1, 1], upper[1, 1]; length = 20),
        range(lower[end, 1], upper[end, 1]; length = 20),
    )
    y = vcat(
        upper[:, 2],
        lower[:, 2],
        range(lower[1, 2], upper[1, 2]; length = 20),
        range(lower[end, 2], upper[end, 2]; length = 20),
    )
    return FoilGeometry(
        upper,
        lower,
        extra_length,
        intra_length,
        front_length,
        back_length,
        perimeter,
        area,
        x,
        y,
    )
end

function _sorted_surface(points::AbstractVector{<:AbstractVector{<:Real}})
    sorted = sort(points; by = first)
    xs = Float64[]
    ys = Float64[]
    i = 1
    while i <= length(sorted)
        j = i
        while j < length(sorted) && sorted[j+1][1] == sorted[i][1]
            j += 1
        end
        push!(xs, Float64(sorted[i][1]))
        push!(ys, sum(point[2] for point in sorted[i:j]) / (j - i + 1))
        i = j + 1
    end
    return hcat(xs, ys)
end

function _offset_surface!(surface::Matrix{Float64}, distance::Float64, direction::Int)
    points = copy(surface)
    for i in axes(surface, 1)
        before = points[max(i - 1, 1), :]
        after = points[min(i + 1, size(points, 1)), :]
        tangent = after - before
        magnitude = hypot(tangent...)
        magnitude > 0 || continue
        normal = [-tangent[2], tangent[1]] ./ magnitude
        surface[i, :] .+= direction .* distance .* normal
    end
    return surface
end
