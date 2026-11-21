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
    length(x) == length(y) || throw(DimensionMismatch("x and y must have the same length"))
    all(isfinite, x) && all(isfinite, y) ||
        throw(ArgumentError("x and y must contain only finite values"))
    length(x) < 2 && return 0.0
    return sum(
        (Float64(y[i]) + Float64(y[i+1])) * (Float64(x[i+1]) - Float64(x[i])) / 2 for
        i = 1:(length(x)-1)
    )
end

"""
    trapz(y) -> Float64

Trapezoidal integration with unit spacing.
"""
function trapz(y::AbstractVector)
    return trapz(0:(length(y)-1), y)
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
    length(x) == length(y) || throw(DimensionMismatch("x and y must have the same length"))
    all(isfinite, x) && all(isfinite, y) ||
        throw(ArgumentError("x and y must contain only finite values"))
    result = zeros(Float64, length(y))
    for i = 2:length(y)
        result[i] =
            result[i-1] +
            (Float64(y[i-1]) + Float64(y[i])) * (Float64(x[i]) - Float64(x[i-1])) / 2
    end
    return result
end

"""
    cumtrapz(y) -> Vector{Float64}

Cumulative trapezoidal integration with unit spacing.
"""
function cumtrapz(y::AbstractVector)
    return cumtrapz(0:(length(y)-1), y)
end

# ---------------------------------------------------------------------------
# Quadratic root (replaces Matlab's roots([a b c]))
# ---------------------------------------------------------------------------

"""
    quadratic_roots(a, b, c) -> (r1, r2)

Return the two roots of a·x² + b·x + c = 0.
"""
function quadratic_roots(a::Real, b::Real, c::Real)
    all(isfinite, (a, b, c)) ||
        throw(ArgumentError("quadratic coefficients must be finite"))
    a != 0 || throw(DomainError(a, "the leading coefficient must be nonzero"))

    discriminant = ComplexF64(b)^2 - 4 * ComplexF64(a) * ComplexF64(c)
    root = sqrt(discriminant)
    q = -0.5 * (ComplexF64(b) + (real(ComplexF64(b)) >= 0 ? root : -root))
    if iszero(q)
        repeated = -ComplexF64(b) / (2 * ComplexF64(a))
        return (repeated, repeated)
    end
    return (q / ComplexF64(a), ComplexF64(c) / q)
end
