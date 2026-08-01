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
    error("Not implemented")
end

"""
    area_inertia_tensor(thick, pts) -> Matrix{Float64}

Second moment of area tensor for a thin element of constant thickness `thick`
traced along the polyline `pts` (N×2 matrix [x y]).

Returned as a 2×2 symmetric matrix [Ix Ixy; Ixy Iy].
"""
function area_inertia_tensor(thick::Float64, pts::Matrix{Float64})
    error("Not implemented")
end

"""
    first_moment_of_area(pts, thick) -> (Qx, Qy)

First moments of area for a thin element of thickness `thick` along `pts`.
"""
function first_moment_of_area(pts::Matrix{Float64}, thick::Float64)
    error("Not implemented")
end

"""
    section_centroid(pts) -> Vector{Float64}

Arc-length centroid of a polyline `pts` (N×2).
"""
function section_centroid(pts::Matrix{Float64})
    error("Not implemented")
end
