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
    A  = sum(areas)
    Sx = sum(areas .* positions[:, 1])
    Sy = sum(areas .* positions[:, 2])
    return [Sx, Sy] ./ A
end

"""
    area_inertia_tensor(thick, pts) -> Matrix{Float64}

Second moment of area tensor for a thin element of constant thickness `thick`
traced along the polyline `pts` (N×2 matrix [x y]).

Returned as a 2×2 symmetric matrix [Ix Ixy; Ixy Iy].
"""
function area_inertia_tensor(thick::Float64, pts::Matrix{Float64})
    dx = diff(pts[:, 1])
    dy = diff(pts[:, 2])
    dl = sqrt.(dx.^2 .+ dy.^2)
    mid_x = pts[1:end-1, 1] .+ dx ./ 2
    mid_y = pts[1:end-1, 2] .+ dy ./ 2
    Ix  = trapz(mid_y.^2 .* dl)
    Iy  = trapz(mid_x.^2 .* dl)
    Ixy = trapz(mid_x .* mid_y .* dl)
    return [Ix Ixy; Ixy Iy] .* thick
end

"""
    first_moment_of_area(pts, thick) -> (Qx, Qy)

First moments of area for a thin element of thickness `thick` along `pts`.
"""
function first_moment_of_area(pts::Matrix{Float64}, thick::Float64)
    dx = diff(pts[:, 1])
    dy = diff(pts[:, 2])
    dl = sqrt.(dx.^2 .+ dy.^2)
    dA = dl .* thick
    mid_x = pts[1:end-1, 1] .+ dx ./ 2
    mid_y = pts[1:end-1, 2] .+ dy ./ 2
    Qx = trapz(mid_x .* dA)
    Qy = trapz(mid_y .* dA)
    return Qx, Qy
end

"""
    section_centroid(pts) -> Vector{Float64}

Arc-length centroid of a polyline `pts` (N×2).
"""
function section_centroid(pts::Matrix{Float64})
    dx = diff(pts[:, 1])
    dy = diff(pts[:, 2])
    dl = sqrt.(dx.^2 .+ dy.^2)
    mid_x = pts[1:end-1, 1] .+ dx ./ 2
    mid_y = pts[1:end-1, 2] .+ dy ./ 2
    L   = sum(dl)
    x_m = trapz(mid_x .* dl) / L
    y_m = trapz(mid_y .* dl) / L
    return [x_m, y_m]
end
