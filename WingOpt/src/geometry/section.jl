"""
Cross-section geometry helpers.
Translates `geometry/section/`: multielement_centroid.m,
area_inertia_tensor.m, first_moment_of_area.m, section_centroid.m.
"""

"""
    multielement_centroid(positions, areas) -> Vector{Float64}

Centroid of a multi-element section by the decomposition method.

`positions`: N×2 matrix of element centroid coordinates [x y].
`areas`:     N-vector of element areas.
"""
function multielement_centroid(positions::Matrix{Float64}, areas::Vector{Float64})
    size(positions, 2) == 2 && size(positions, 1) == length(areas) ||
        throw(DimensionMismatch("positions must be N×2 and areas must have length N"))
    total_area = sum(areas)
    total_area > 0 || throw(ArgumentError("total area must be positive"))
    return vec(sum(positions .* areas; dims = 1) ./ total_area)
end

"""
    area_inertia_tensor(thick, pts) -> Matrix{Float64}

Second moment of area tensor for a thin element of constant thickness `thick`
traced along the polyline `pts` (N×2 matrix [x y]).

Returned as a 2×2 symmetric matrix [Ix Ixy; Ixy Iy].
"""
function area_inertia_tensor(thick::Float64, pts::Matrix{Float64})
    size(pts, 2) == 2 || throw(DimensionMismatch("points must be N×2"))
    size(pts, 1) >= 2 || return zeros(2, 2)
    x0, y0 = pts[1:(end-1), 1], pts[1:(end-1), 2]
    x1, y1 = pts[2:end, 1], pts[2:end, 2]
    dx, dy = x1 .- x0, y1 .- y0
    dl = hypot.(dx, dy)
    Ix = thick * sum(dl .* (y0 .^ 2 .+ y0 .* y1 .+ y1 .^ 2) ./ 3)
    Iy = thick * sum(dl .* (x0 .^ 2 .+ x0 .* x1 .+ x1 .^ 2) ./ 3)
    Ixy = thick * sum(dl .* (2x0 .* y0 .+ x0 .* y1 .+ x1 .* y0 .+ 2x1 .* y1) ./ 6)
    return [Ix Ixy; Ixy Iy]
end

"""
    first_moment_of_area(pts, thick) -> (Qx, Qy)

First moments of area for a thin element of thickness `thick` along `pts`.
"""
function first_moment_of_area(pts::Matrix{Float64}, thick::Float64)
    size(pts, 2) == 2 || throw(DimensionMismatch("points must be N×2"))
    size(pts, 1) >= 2 || return (0.0, 0.0)
    x0, y0 = pts[1:(end-1), 1], pts[1:(end-1), 2]
    x1, y1 = pts[2:end, 1], pts[2:end, 2]
    dx, dy = x1 .- x0, y1 .- y0
    dl = hypot.(dx, dy)
    return (thick * sum((x0 .+ x1) .* dl ./ 2), thick * sum((y0 .+ y1) .* dl ./ 2))
end

"""
    section_centroid(pts) -> Vector{Float64}

Arc-length centroid of a polyline `pts` (N×2).
"""
function section_centroid(pts::Matrix{Float64})
    size(pts, 2) == 2 || throw(DimensionMismatch("points must be N×2"))
    size(pts, 1) >= 2 || throw(ArgumentError("at least two points are required"))
    dx, dy = diff(pts[:, 1]), diff(pts[:, 2])
    dl = hypot.(dx, dy)
    length_total = sum(dl)
    length_total > 0 || throw(ArgumentError("polyline must have nonzero length"))
    x, y = pts[1:(end-1), 1] .+ dx ./ 2, pts[1:(end-1), 2] .+ dy ./ 2
    return [sum(x .* dl) / length_total, sum(y .* dl) / length_total]
end
