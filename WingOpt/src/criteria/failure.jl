"""
Failure criteria for composite/orthotropic materials.
Translates `criterias/failure/`: max_tension_criteria.m,
tsai_wu_criteria.m, tsai_hill_criteria.m.

All functions accept a 3×3 stress tensor `sigma` and a `MaterialProperties`.
"""

# ---------------------------------------------------------------------------
# Helper: extract strengths based on orientation flag
# ---------------------------------------------------------------------------

function _strengths(mat::MaterialProperties, opt::String)
    if opt == "axial"
        return mat.axial_trac, mat.axial_comp,
               mat.trans_trac, mat.trans_comp, mat.shear
    elseif opt == "trans"
        # swap axial ↔ transverse
        return mat.trans_trac, mat.trans_comp,
               mat.axial_trac, mat.axial_comp, mat.shear
    else
        error("Unknown orientation option: \"$opt\". Use \"axial\" or \"trans\".")
    end
end

# ---------------------------------------------------------------------------
# 1. Maximum-stress criterion
# ---------------------------------------------------------------------------

"""
    max_tension_criteria(sigma, mat, opt) -> Float64

Maximum-stress failure criterion. Returns the minimum margin of safety
(positive = safe, negative = failed).

`sigma` – 3×3 stress tensor  
`mat`   – MaterialProperties  
`opt`   – \"axial\" or \"trans\" (fibre direction)
"""
function max_tension_criteria(sigma::Matrix{Float64},
                               mat::MaterialProperties, opt::String)
    F1t, F1c, F2t, F2c, F6 = _strengths(mat, opt)

    σ1 = sigma[1, 1]
    σ2 = sigma[2, 2]
    τ6 = sigma[1, 2]

    ms = Float64[]
    if σ1 >= 0
        push!(ms, F1t / σ1  - 1)
    elseif σ1 < 0
        push!(ms, -F1c / σ1 - 1)
    end
    if σ2 >= 0
        push!(ms, F2t / σ2  - 1)
    elseif σ2 < 0
        push!(ms, -F2c / σ2 - 1)
    end
    if abs(τ6) > eps()
        push!(ms, abs(F6 / τ6) - 1)
    end

    # Filter ±Inf entries (zero-stress components)
    finite_ms = filter(isfinite, ms)
    isempty(finite_ms) && return Inf
    return minimum(finite_ms)
end

# ---------------------------------------------------------------------------
# 2. Tsai-Wu criterion
# ---------------------------------------------------------------------------

"""
    tsai_wu_criteria(sigma, mat, opt) -> Float64

Tsai-Wu polynomial failure criterion.
Returns the margin of safety: sf - 1, where sf is the safety factor
at failure (positive = safe, negative = failed).
"""
function tsai_wu_criteria(sigma::Matrix{Float64},
                           mat::MaterialProperties, opt::String)
    F1t, F1c, F2t, F2c, F6 = _strengths(mat, opt)

    σ1 = sigma[1, 1];  σ2 = sigma[2, 2];  τ6 = sigma[1, 2]

    f1  = 1/F1t - 1/F1c
    f11 = 1/(F1t * F1c)
    f2  = 1/F2t - 1/F2c
    f22 = 1/(F2t * F2c)
    f66 = 1/F6^2
    f12 = -0.5 * sqrt(f11 * f22)

    a = f11*σ1^2 + f22*σ2^2 + f66*τ6^2 + 2*f12*σ1*σ2
    b = f1*σ1    + f2*σ2

    # Solve  a·sf² + b·sf - 1 = 0
    a < eps() && return (b > 0 ? 1/b - 1 : Inf)   # degenerate: linear
    r1, r2 = quadratic_roots(a, b, -1.0)
    sf = max(r1, r2)
    return sf - 1
end

# ---------------------------------------------------------------------------
# 3. Tsai-Hill criterion
# ---------------------------------------------------------------------------

"""
    tsai_hill_criteria(sigma, mat) -> Float64

Tsai-Hill failure criterion. Returns 1 - FI (positive = safe, 0 = onset of failure).
"""
function tsai_hill_criteria(sigma::Matrix{Float64}, mat::MaterialProperties)
    σ1 = sigma[1, 1];  σ2 = sigma[2, 2];  τ6 = sigma[1, 2]

    F1 = σ1 >= 0 ? mat.axial_trac : mat.axial_comp
    F2 = σ2 >= 0 ? mat.trans_trac : mat.trans_comp
    F6 = mat.shear

    FI = (σ1/F1)^2 + (σ2/F2)^2 - σ1*σ2/F1^2 + (τ6/F6)^2
    return 1.0 - FI    # positive = safe (different sign convention from MS)
end
