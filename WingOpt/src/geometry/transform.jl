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
    error("Not implemented")
end

"""
    rotate_inertia_tensor(I, theta) -> Matrix{Float64}

Rotate 2×2 inertia tensor `I` by angle `theta` (radians).
"""
function rotate_inertia_tensor(I::Matrix{Float64}, theta::Float64)
    error("Not implemented")
end

"""
    principal_inertia_tensor(I) -> (Ip, theta_p)

Compute the principal inertia tensor `Ip` and the principal-axis angle `theta_p`
(radians) using Mohr's circle.

Returns `(Ip::Matrix{Float64}, theta_p::Float64)`.
"""
function principal_inertia_tensor(I::Matrix{Float64})
    error("Not implemented")
end
