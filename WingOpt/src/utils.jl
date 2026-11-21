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
    return trapz(collect(0:(length(y)-1)), y)
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
    result = zeros(Float64, length(x))
    for i = 2:length(x)
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
    return cumtrapz(collect(0:(length(y)-1)), y)
end

# ---------------------------------------------------------------------------
# Quadratic root (replaces Matlab's roots([a b c]))
# ---------------------------------------------------------------------------

"""
    quadratic_roots(a, b, c) -> (r1, r2)

Return the two roots of a·x² + b·x + c = 0.
"""
function quadratic_roots(a::Float64, b::Float64, c::Float64)
    error("Not implemented")
end
