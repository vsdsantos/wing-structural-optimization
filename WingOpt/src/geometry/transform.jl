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
    size(I) == (2, 2) && length(d) == 2 || throw(DimensionMismatch("expected a 2×2 tensor and a 2-vector"))
    return [I[1, 1] + d[1]^2 * A I[1, 2] + d[1] * d[2] * A;
            I[2, 1] + d[1] * d[2] * A I[2, 2] + d[2]^2 * A]
end

"""
    rotate_inertia_tensor(I, theta) -> Matrix{Float64}

Rotate 2×2 inertia tensor `I` by angle `theta` (radians).
"""
function rotate_inertia_tensor(I::Matrix{Float64}, theta::Float64)
    size(I) == (2, 2) || throw(DimensionMismatch("inertia tensor must be 2×2"))
    c, s = cos(2theta), sin(2theta)
    meanI = (I[1, 1] + I[2, 2]) / 2
    return [meanI + (I[2, 2] - I[1, 1]) * c / 2 + I[1, 2] * s
            (I[1, 1] - I[2, 2]) * s / 2 + I[1, 2] * c;
            (I[1, 1] - I[2, 2]) * s / 2 + I[2, 1] * c
            meanI + (I[1, 1] - I[2, 2]) * c / 2 - I[2, 1] * s]
end

"""
    principal_inertia_tensor(I) -> (Ip, theta_p)

Compute the principal inertia tensor `Ip` and the principal-axis angle `theta_p`
(radians) using Mohr's circle.

Returns `(Ip::Matrix{Float64}, theta_p::Float64)`.
"""
function principal_inertia_tensor(I::Matrix{Float64})
    size(I) == (2, 2) || throw(DimensionMismatch("inertia tensor must be 2×2"))
    I[1, 2] ≈ I[2, 1] || throw(ArgumentError("inertia tensor must be symmetric"))
    meanI = (I[1, 1] + I[2, 2]) / 2
    radius = hypot((I[1, 1] - I[2, 2]) / 2, I[1, 2])
    angle = iszero(radius) ? 0.0 : atan(2I[1, 2], I[1, 1] - I[2, 2]) / 2
    return ([meanI + radius 0.0; 0.0 meanI - radius], angle)
end
