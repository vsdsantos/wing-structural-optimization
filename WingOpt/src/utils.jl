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
    length(x) == 1 && return 0.0
    result = 0.0
    @inbounds for i in 2:length(x)
        result += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    return result / 2
end

"""
    trapz(y) -> Float64

Trapezoidal integration with unit spacing.
"""
trapz(y::AbstractVector) = trapz(eachindex(y), y)

# ---------------------------------------------------------------------------
# Cumulative trapezoidal integration
# ---------------------------------------------------------------------------

"""
    cumtrapz(x, y) -> Vector{Float64}

Cumulative trapezoidal integration of `y` over `x`.
Equivalent to Matlab's `cumtrapz(x, y)`.
"""
function cumtrapz(x::AbstractVector, y::AbstractVector)
    n = length(x)
    result = zeros(n)
    @inbounds for i in 2:n
        result[i] = result[i-1] + (x[i] - x[i-1]) * (y[i] + y[i-1]) / 2
    end
    return result
end

"""
    cumtrapz(y) -> Vector{Float64}

Cumulative trapezoidal integration with unit spacing.
"""
cumtrapz(y::AbstractVector) = cumtrapz(eachindex(y), y)

# ---------------------------------------------------------------------------
# Quadratic root (replaces Matlab's roots([a b c]))
# ---------------------------------------------------------------------------

"""
    quadratic_roots(a, b, c) -> (r1, r2)

Return the two roots of a·x² + b·x + c = 0.
"""
function quadratic_roots(a::Float64, b::Float64, c::Float64)
    disc = b^2 - 4*a*c
    sq = sqrt(max(disc, 0.0))
    return (-b + sq) / (2a), (-b - sq) / (2a)
end
