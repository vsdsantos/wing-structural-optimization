"""
Inertia-tensor operations.
Translates `geometry/transform/`: translate_inertia_tensor.m,
rotate_inertia_tensor.m, principal_inertia_tensor.m.
"""

"""
    translate_inertia_tensor(I, d, A) -> Matrix{Float64}

Parallel-axis theorem: shift a 2×2 inertia tensor `I` (m⁴) by displacement
`d = [dx, dy]` (m) for an element of area `A` (m²).
Equivalent to Matlab's `translada_eixo` / `translate_inertia_tensor`.
"""
function translate_inertia_tensor(I::Matrix{Float64}, d::Vector{Float64}, A::Float64)
    _validate_inertia_tensor(I)
    length(d) == 2 || throw(DimensionMismatch("displacement must contain [dx, dy]"))
    all(isfinite, I) && all(isfinite, d) && isfinite(A) ||
        throw(ArgumentError("inertia tensor, displacement, and area must be finite"))
    A >= 0 || throw(DomainError(A, "area must be nonnegative"))
    dx, dy = d
    return [
        I[1, 1] + dx^2 * A I[1, 2] + dx * dy * A;
        I[2, 1] + dx * dy * A I[2, 2] + dy^2 * A
    ]
end

"""
    rotate_inertia_tensor(I, theta) -> Matrix{Float64}

Rotate 2×2 inertia tensor `I` (m⁴) by angle `theta` (radians).
"""
function rotate_inertia_tensor(I::Matrix{Float64}, theta::Float64)
    _validate_inertia_tensor(I)
    all(isfinite, I) && isfinite(theta) ||
        throw(ArgumentError("inertia tensor and angle must be finite"))
    c, s = cos(2theta), sin(2theta)
    mean_inertia = (I[1, 1] + I[2, 2]) / 2
    difference = (I[1, 1] - I[2, 2]) / 2
    product = (I[1, 2] + I[2, 1]) / 2
    return [
        mean_inertia - difference * c + product * s difference * s + product * c;
        difference * s + product * c mean_inertia + difference * c - product * s
    ]
end

"""
    principal_inertia_tensor(I) -> (Ip, theta_p)

Compute the principal inertia tensor `Ip` and the principal-axis angle `theta_p`
(radians) using Mohr's circle.

Returns `(Ip::Matrix{Float64}, theta_p::Float64)`, with the larger principal
moment in `Ip[1, 1]`; rotating `I` by `theta_p` gives `Ip`.
"""
function principal_inertia_tensor(I::Matrix{Float64})
    _validate_inertia_tensor(I)
    theta = atan(2 * I[1, 2], I[2, 2] - I[1, 1]) / 2
    rotated = rotate_inertia_tensor(I, theta)
    principal = [rotated[1, 1] 0.0; 0.0 rotated[2, 2]]
    return principal, theta
end

function _validate_inertia_tensor(I::Matrix{Float64})
    size(I) == (2, 2) || throw(DimensionMismatch("inertia tensor must be 2×2"))
    all(isfinite, I) || throw(ArgumentError("inertia tensor must be finite"))
    isapprox(I, I'; rtol = 1e-10, atol = 1e-14) ||
        throw(ArgumentError("inertia tensor must be symmetric"))
    minimum(eigvals(Symmetric(I))) >= -1e-12 ||
        throw(DomainError(I, "inertia tensor must be positive semidefinite"))
    return nothing
end
