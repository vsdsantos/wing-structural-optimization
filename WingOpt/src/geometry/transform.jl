"""
Inertia-tensor operations.
Translates `geometry/transform/`: translate_inertia_tensor.m,
rotate_inertia_tensor.m, principal_inertia_tensor.m.
"""

"""
    translate_inertia_tensor(I, d, A) -> Matrix{Float64}

Parallel-axis theorem: shift a 2×2 inertia tensor `I` by displacement `d = [dx, dy]`
for an element of area `A`.
Equivalent to Matlab's `translada_eixo` / `translate_inertia_tensor`.
"""
function translate_inertia_tensor(I::Matrix{Float64}, d::Vector{Float64}, A::Float64)
    size(I) == (2, 2) || throw(DimensionMismatch("inertia tensor must be 2×2"))
    length(d) == 2 || throw(DimensionMismatch("displacement must have length 2"))
    dx, dy = d
    return I .+ A .* [dx^2 dx * dy; dx * dy dy^2]
end

"""
    rotate_inertia_tensor(I, theta) -> Matrix{Float64}

Rotate 2×2 inertia tensor `I` by angle `theta` (radians).
"""
function rotate_inertia_tensor(I::Matrix{Float64}, theta::Float64)
    size(I) == (2, 2) || throw(DimensionMismatch("inertia tensor must be 2×2"))
    c, s = cos(2theta), sin(2theta)
    Ixx, Iyy, Ixy, Iyx = I[1, 1], I[2, 2], I[1, 2], I[2, 1]
    rotated_xx = (Ixx + Iyy) / 2 + (Iyy - Ixx) * c / 2 + Ixy * s
    rotated_yy = (Ixx + Iyy) / 2 + (Ixx - Iyy) * c / 2 - Iyx * s
    rotated_xy = (Ixx - Iyy) * s / 2 + Ixy * c
    rotated_yx = (Ixx - Iyy) * s / 2 + Iyx * c
    return [rotated_xx rotated_xy; rotated_yx rotated_yy]
end

"""
    principal_inertia_tensor(I) -> (Ip, theta_p)

Compute the principal inertia tensor `Ip` and the principal-axis angle `theta_p`
(radians) using Mohr's circle.

Returns `(Ip::Matrix{Float64}, theta_p::Float64)`.
"""
function principal_inertia_tensor(I::Matrix{Float64})
    size(I) == (2, 2) || throw(DimensionMismatch("inertia tensor must be 2×2"))
    mean_i = tr(I) / 2
    radius = hypot((I[1, 1] - I[2, 2]) / 2, I[1, 2])
    angle = iszero(I[2, 2] - I[1, 1]) ?
            (iszero(I[1, 2]) ? 0.0 : pi / 4) :
            0.5 * atan(abs(2I[1, 2] / (I[1, 1] - I[2, 2])))
    return ([mean_i + radius 0.0; 0.0 mean_i - radius], angle)
end
