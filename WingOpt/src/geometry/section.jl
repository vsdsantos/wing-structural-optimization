"""
Cross-section geometry helpers.
Translates `geometry/section/`: multielement_centroid.m,
area_inertia_tensor.m, first_moment_of_area.m, section_centroid.m.
"""

"""
    multielement_centroid(positions, areas) -> Vector{Float64}

Centroid of a multi-element section by the decomposition method.

`positions`: N×2 matrix of element centroid coordinates [x y] in metres.
`areas`:     N-vector of element areas in m². Returns the centroid in metres.
"""
function multielement_centroid(positions::Matrix{Float64}, areas::Vector{Float64})
    size(positions, 2) == 2 ||
        throw(DimensionMismatch("positions must have two columns [x, y]"))
    size(positions, 1) == length(areas) ||
        throw(DimensionMismatch("one area is required for each position"))
    all(isfinite, positions) && all(isfinite, areas) ||
        throw(ArgumentError("positions and areas must be finite"))
    all(>=(0), areas) || throw(DomainError(areas, "areas must be nonnegative"))
    total_area = sum(areas)
    total_area > 0 || throw(DomainError(total_area, "total area must be positive"))
    return vec(sum(positions .* reshape(areas, :, 1); dims = 1) ./ total_area)
end

"""
    area_inertia_tensor(thick, pts) -> Matrix{Float64}

Second moment of area tensor (m⁴) for a thin element of thickness `thick` (m)
traced along the polyline `pts` (N×2 matrix [x y] in metres).

Returned as a 2×2 symmetric matrix [Ix Ixy; Ixy Iy].
"""
function area_inertia_tensor(thick::Float64, pts::Matrix{Float64})
    _validate_polyline(pts)
    isfinite(thick) && thick >= 0 ||
        throw(DomainError(thick, "thickness must be finite and nonnegative"))
    x0, x1 = pts[1:(end-1), 1], pts[2:end, 1]
    y0, y1 = pts[1:(end-1), 2], pts[2:end, 2]
    dx, dy = x1 - x0, y1 - y0
    dl = hypot.(dx, dy)
    Ix = thick * sum((y0 .^ 2 .+ y0 .* y1 .+ y1 .^ 2) .* dl ./ 3)
    Iy = thick * sum((x0 .^ 2 .+ x0 .* x1 .+ x1 .^ 2) .* dl ./ 3)
    Ixy = thick * sum((2 .* x0 .* y0 .+ x0 .* y1 .+ x1 .* y0 .+ 2 .* x1 .* y1) .* dl ./ 6)
    return [Ix Ixy; Ixy Iy]
end

"""
    first_moment_of_area(pts, thick) -> (Qx, Qy)

First moments of area (m³) for a thin element of thickness `thick` (m) along
`pts` (coordinates in metres).
"""
function first_moment_of_area(pts::Matrix{Float64}, thick::Float64)
    _validate_polyline(pts)
    isfinite(thick) && thick >= 0 ||
        throw(DomainError(thick, "thickness must be finite and nonnegative"))
    x0, x1 = pts[1:(end-1), 1], pts[2:end, 1]
    y0, y1 = pts[1:(end-1), 2], pts[2:end, 2]
    dx, dy = x1 - x0, y1 - y0
    dl = hypot.(dx, dy)
    Qx = thick * sum((y0 .+ y1) .* dl ./ 2)
    Qy = thick * sum((x0 .+ x1) .* dl ./ 2)
    return Qx, Qy
end

"""
    section_centroid(pts) -> Vector{Float64}

Arc-length centroid (metres) of a polyline `pts` (N×2, coordinates in metres).
"""
function section_centroid(pts::Matrix{Float64})
    _validate_polyline(pts)
    x0, x1 = pts[1:(end-1), 1], pts[2:end, 1]
    y0, y1 = pts[1:(end-1), 2], pts[2:end, 2]
    dx, dy = x1 - x0, y1 - y0
    dl = hypot.(dx, dy)
    total_length = sum(dl)
    total_length > 0 || throw(DomainError(total_length, "polyline length must be positive"))
    return [sum((x0 .+ x1) .* dl ./ 2), sum((y0 .+ y1) .* dl ./ 2)] ./ total_length
end

function _validate_polyline(pts::Matrix{Float64})
    size(pts, 2) == 2 || throw(DimensionMismatch("polyline must have two columns [x, y]"))
    size(pts, 1) >= 2 || throw(ArgumentError("polyline must contain at least two points"))
    all(isfinite, pts) || throw(ArgumentError("polyline coordinates must be finite"))
    return nothing
end
