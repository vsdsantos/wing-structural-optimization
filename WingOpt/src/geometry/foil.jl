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
    path = joinpath(
        foil_dir,
        endswith(lowercase(foil_name), ".dat") ? foil_name : "$foil_name.dat",
    )
    isfile(path) || throw(ArgumentError("airfoil file not found: $path"))
    points = NTuple{2,Float64}[]
    for (line_number, line) in enumerate(eachline(path))
        fields = split(strip(line))
        length(fields) < 2 && continue
        x, y = tryparse(Float64, fields[1]), tryparse(Float64, fields[2])
        (x === nothing || y === nothing) && continue
        all(isfinite, (x, y)) ||
            throw(ArgumentError("non-finite coordinate on line $line_number of $path"))
        push!(points, (x, y))
    end
    length(points) >= 4 ||
        throw(ArgumentError("airfoil file must contain at least four coordinate points"))
    data = Matrix{Float64}(undef, length(points), 2)
    for (i, point) in enumerate(points)
        data[i, :] .= point
    end
    leading_edge = argmin(data[:, 1])
    upper = data[1:leading_edge, :]
    lower = data[leading_edge:end, :]
    upper = upper[sortperm(upper[:, 1]), :]
    lower = lower[sortperm(lower[:, 1]), :]
    function unique_surface(surface)
        xs, ys = Float64[], Float64[]
        for i in eachindex(surface[:, 1])
            if !isempty(xs) && isapprox(surface[i, 1], xs[end]; atol = 1e-12)
                ys[end] = (ys[end] + surface[i, 2]) / 2
            else
                push!(xs, surface[i, 1])
                push!(ys, surface[i, 2])
            end
        end
        length(xs) >= 2 || throw(
            ArgumentError("each airfoil surface needs at least two distinct x coordinates"),
        )
        return xs, ys
    end
    ux, uy = unique_surface(upper)
    lx, ly = unique_surface(lower)
    x_min, x_max = minimum(data[:, 1]), maximum(data[:, 1])
    ux .= (ux .- x_min) ./ (x_max - x_min)
    lx .= (lx .- x_min) ./ (x_max - x_min)
    minimum(diff(ux)) > 0 && minimum(diff(lx)) > 0 ||
        throw(ArgumentError("airfoil x coordinates must increase"))
    upper_interp = _linear_interpolant(ux, uy)
    lower_interp = _linear_interpolant(lx, ly)
    spacing = min(minimum(diff(ux)), minimum(diff(lx))) * 0.1
    return FoilData(upper_interp, lower_interp, hcat(ux, uy), hcat(lx, ly), spacing)
end

function _linear_interpolant(x::Vector{Float64}, y::Vector{Float64})
    return function (value)
        x0 = Float64(value)
        x[1] - 1e-12 <= x0 <= x[end] + 1e-12 || throw(BoundsError(x, x0))
        x0 = clamp(x0, x[1], x[end])
        i = clamp(searchsortedlast(x, x0), 1, length(x) - 1)
        return y[i] + (x0 - x[i]) * (y[i+1] - y[i]) / (x[i+1] - x[i])
    end
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
    chord > 0 || throw(ArgumentError("chord must be positive"))
    length(limits) == 2 && 0 <= limits[1] < limits[2] <= 1 ||
        throw(ArgumentError("limits must satisfy 0 ≤ front < back ≤ 1"))
    length(offset) == 2 || throw(DimensionMismatch("offset must have two values"))
    step = foil.esp > 0 ? foil.esp : 0.01
    xs = collect(range(Float64(limits[1]), Float64(limits[2]); step = step))
    isempty(xs) && (xs = [Float64(limits[1])])
    xs[end] < limits[2] - eps(Float64) && push!(xs, Float64(limits[2]))
    x_full = range(0.0, 1.0; length = max(101, length(xs)))
    camber = (foil.ex_lin(0.25) + foil.in_lin(0.25)) / 2
    shift = camber * chord
    upper_y = foil.ex_lin.(xs)
    lower_y = foil.in_lin.(xs)
    upper = hcat((xs .- 0.25) .* chord, upper_y .* chord .- shift .- offset[1])
    lower = hcat((xs .- 0.25) .* chord, lower_y .* chord .- shift .+ offset[2])
    if !iszero(scale_thick)
        _offset_surface!(upper, scale_thick)
        _offset_surface!(lower, -scale_thick)
    end
    segment_length(points) = sum(hypot.(diff(points[:, 1]), diff(points[:, 2])))
    front_length = upper[1, 2] - lower[1, 2]
    back_length = upper[end, 2] - lower[end, 2]
    area = trapz(upper[:, 1], upper[:, 2] .- lower[:, 2])
    x_closed = vcat(
        upper[:, 1],
        reverse(lower[:, 1]),
        fill(upper[1, 1], 2),
        fill(upper[end, 1], 2),
    )
    y_closed = vcat(
        upper[:, 2],
        reverse(lower[:, 2]),
        [lower[1, 2], upper[1, 2], upper[end, 2], lower[end, 2]],
    )
    return FoilGeometry(
        upper,
        lower,
        segment_length(upper),
        segment_length(lower),
        front_length,
        back_length,
        segment_length(upper) +
        segment_length(lower) +
        abs(front_length) +
        abs(back_length),
        area,
        x_closed,
        y_closed,
    )
end

function _offset_surface!(points::Matrix{Float64}, offset::Float64)
    for i in axes(points, 1)
        lo, hi = max(1, i - 1), min(size(points, 1), i + 1)
        tangent = points[hi, :] .- points[lo, :]
        normal = [-tangent[2], tangent[1]]
        norm(normal) > 0 && (points[i, :] .+= offset .* normal ./ norm(normal))
    end
    return points
end
