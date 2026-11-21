"""
Numerical utility functions.
Provides trapz / cumtrapz equivalents (Matlab-compatible signatures) without
requiring an external integration package.
"""

# ---------------------------------------------------------------------------
# Trapezoidal integration
# ---------------------------------------------------------------------------

"""
    trapz(x, y) -> Float64

Trapezoidal integration of `y` over `x`.
Equivalent to Matlab's `trapz(x, y)`.
"""
function trapz(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || throw(DimensionMismatch("x and y must have equal lengths"))
    length(x) < 2 && return 0.0
    return sum((x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2 for i in 1:(length(x) - 1))
end

"""
    trapz(y) -> Float64

Trapezoidal integration with unit spacing.
"""
function trapz(y::AbstractVector)
    return trapz(1:length(y), y)
end

# ---------------------------------------------------------------------------
# Cumulative trapezoidal integration
# ---------------------------------------------------------------------------

"""
    cumtrapz(x, y) -> Vector{Float64}

Cumulative trapezoidal integration of `y` over `x`.
Equivalent to Matlab's `cumtrapz(x, y)`.
"""
function cumtrapz(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || throw(DimensionMismatch("x and y must have equal lengths"))
    result = zeros(Float64, length(y))
    for i in 2:length(y)
        result[i] = result[i - 1] + (x[i] - x[i - 1]) * (y[i] + y[i - 1]) / 2
    end
    return result
end

"""
    cumtrapz(y) -> Vector{Float64}

Cumulative trapezoidal integration with unit spacing.
"""
function cumtrapz(y::AbstractVector)
    return cumtrapz(1:length(y), y)
end

# ---------------------------------------------------------------------------
# Quadratic root (replaces Matlab's roots([a b c]))
# ---------------------------------------------------------------------------

"""
    quadratic_roots(a, b, c) -> (r1, r2)

Return the two roots of a·x² + b·x + c = 0.
"""
function quadratic_roots(a::Float64, b::Float64, c::Float64)
    if iszero(a)
        iszero(b) && throw(ArgumentError("the equation has no unique roots"))
        root = -c / b
        return (root, root)
    end
    discriminant = b^2 - 4a * c
    discriminant < 0 && return (complex(-b, sqrt(-discriminant)) / (2a),
                                complex(-b, -sqrt(-discriminant)) / (2a))
    root = sqrt(discriminant)
    q = -0.5 * (b + copysign(root, b))
    iszero(q) && return (0.0, 0.0)
    return (q / a, c / q)
end
