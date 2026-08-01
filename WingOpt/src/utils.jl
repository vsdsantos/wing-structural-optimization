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
    error("Not implemented")
end

"""
    trapz(y) -> Float64

Trapezoidal integration with unit spacing.
"""
function trapz(y::AbstractVector)
    error("Not implemented")
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
    error("Not implemented")
end

"""
    cumtrapz(y) -> Vector{Float64}

Cumulative trapezoidal integration with unit spacing.
"""
function cumtrapz(y::AbstractVector)
    error("Not implemented")
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
