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
    size(positions, 2) == 2 || throw(DimensionMismatch("positions must have two columns"))
    size(positions, 1) == length(areas) ||
        throw(DimensionMismatch("one area is required per position"))
    total_area = sum(areas)
    iszero(total_area) && throw(DomainError(total_area, "total area must be nonzero"))
    return vec(sum(positions .* areas; dims = 1) ./ total_area)
end

"""
    area_inertia_tensor(thick, pts) -> Matrix{Float64}

Second moment of area tensor for a thin element of constant thickness `thick`
traced along the polyline `pts` (N×2 matrix [x y]).

Returned as a 2×2 symmetric matrix [Ix Ixy; Ixy Iy].
"""
function area_inertia_tensor(thick::Float64, pts::Matrix{Float64})
    size(pts, 2) == 2 || throw(DimensionMismatch("points must have two columns"))
    thick >= 0 || throw(DomainError(thick, "thickness must be nonnegative"))
    size(pts, 1) < 2 && return zeros(2, 2)
    delta = diff(pts; dims = 1)
    segment_length = hypot.(delta[:, 1], delta[:, 2])
    midpoint = (pts[1:(end - 1), :] .+ pts[2:end, :]) ./ 2
    return thick .* [
        sum(midpoint[:, 2] .^ 2 .* segment_length) sum(
            midpoint[:, 1] .* midpoint[:, 2] .* segment_length
        )
        sum(midpoint[:, 1] .* midpoint[:, 2] .* segment_length) sum(
            midpoint[:, 1] .^ 2 .* segment_length
        )
    ]
end

"""
    first_moment_of_area(pts, thick) -> (Qx, Qy)

First moments of area for a thin element of thickness `thick` along `pts`.
"""
function first_moment_of_area(pts::Matrix{Float64}, thick::Float64)
    size(pts, 2) == 2 || throw(DimensionMismatch("points must have two columns"))
    thick >= 0 || throw(DomainError(thick, "thickness must be nonnegative"))
    size(pts, 1) < 2 && return (0.0, 0.0)
    delta = diff(pts; dims = 1)
    segment_length = hypot.(delta[:, 1], delta[:, 2])
    midpoint = (pts[1:(end - 1), :] .+ pts[2:end, :]) ./ 2
    area = thick .* segment_length
    return (sum(midpoint[:, 1] .* area), sum(midpoint[:, 2] .* area))
end

"""
    section_centroid(pts) -> Vector{Float64}

Arc-length centroid of a polyline `pts` (N×2).
"""
function section_centroid(pts::Matrix{Float64})
    size(pts, 2) == 2 || throw(DimensionMismatch("points must have two columns"))
    size(pts, 1) < 2 && throw(ArgumentError("at least two points are required"))
    delta = diff(pts; dims = 1)
    segment_length = hypot.(delta[:, 1], delta[:, 2])
    total_length = sum(segment_length)
    iszero(total_length) &&
        throw(DomainError(total_length, "polyline length must be nonzero"))
    midpoint = (pts[1:(end - 1), :] .+ pts[2:end, :]) ./ 2
    return vec(sum(midpoint .* segment_length; dims = 1) ./ total_length)
end
