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
    return [I[1,1] + d[1]^2 * A    prod(d)*A + I[1,2];
            prod(d)*A + I[2,1]    I[2,2] + d[2]^2 * A]
end

"""
    rotate_inertia_tensor(I, theta) -> Matrix{Float64}

Rotate 2×2 inertia tensor `I` by angle `theta` (radians).
"""
function rotate_inertia_tensor(I::Matrix{Float64}, theta::Float64)
    Ir = zeros(2, 2)
    c2 = cos(2*theta);  s2 = sin(2*theta)
    Ir[1,1] =  0.5*(I[1,1]+I[2,2]) + 0.5*(I[2,2]-I[1,1])*c2 + I[1,2]*s2
    Ir[2,2] =  0.5*(I[1,1]+I[2,2]) + 0.5*(I[1,1]-I[2,2])*c2 - I[2,1]*s2
    Ir[1,2] =  0.5*(I[1,1]-I[2,2])*s2 + I[1,2]*c2
    Ir[2,1] =  0.5*(I[1,1]-I[2,2])*s2 + I[2,1]*c2
    return Ir
end

"""
    principal_inertia_tensor(I) -> (Ip, theta_p)

Compute the principal inertia tensor `Ip` and the principal-axis angle `theta_p`
(radians) using Mohr's circle.

Returns `(Ip::Matrix{Float64}, theta_p::Float64)`.
"""
function principal_inertia_tensor(I::Matrix{Float64})
    Imed   = (I[1,1] + I[2,2]) / 2
    denom  = I[2,2] - I[1,1]
    R      = sqrt((0.5 * denom)^2 + I[1,2]^2)

    Ip = zeros(2, 2)
    Ip[1,1] = Imed + R
    Ip[2,2] = Imed - R

    # guard against division by zero when the section is already principal
    theta_p = abs(denom) < eps(Float64) ? 0.0 :
              0.5 * atan(abs(2 * I[1,2] / denom))

    return Ip, theta_p
end
