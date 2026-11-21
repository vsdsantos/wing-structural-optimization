"""
Failure criteria for composite/orthotropic materials.
Translates `criterias/failure/`: max_tension_criteria.m,
tsai_wu_criteria.m, tsai_hill_criteria.m.

All functions accept a 3×3 stress tensor `sigma` and a `MaterialProperties`.
"""

"""
    max_tension_criteria(sigma, mat, opt) -> Float64

Maximum-stress failure criterion. Returns the minimum margin of safety
(positive = safe, negative = failed).

`sigma` – 3×3 stress tensor  
`mat`   – MaterialProperties  
`opt`   – \"axial\" or \"trans\" (fibre direction)
"""
function max_tension_criteria(sigma::Matrix{Float64}, mat::MaterialProperties, opt::String)
    size(sigma) == (3, 3) || throw(DimensionMismatch("stress tensor must be 3×3"))
    strengths = _failure_strengths(mat, opt)
    σ1, σ2, τ = sigma[1, 1], sigma[2, 2], sigma[1, 2]
    margin(σ, tensile, compressive) =
        iszero(σ) ? Inf : (σ > 0 ? tensile / σ : -compressive / σ) - 1
    return minimum((
        margin(σ1, strengths[1], strengths[2]),
        margin(σ2, strengths[3], strengths[4]),
        iszero(τ) ? Inf : strengths[5] / abs(τ) - 1,
    ))
end

function _failure_strengths(mat::MaterialProperties, opt::String)
    opt in ("axial", "trans") || throw(ArgumentError("opt must be \"axial\" or \"trans\""))
    axial = (mat.axial_trac, mat.axial_comp, mat.trans_trac, mat.trans_comp)
    transverse = (mat.trans_trac, mat.trans_comp, mat.axial_trac, mat.axial_comp)
    return (opt == "axial" ? axial : transverse)..., mat.shear
end

"""
    tsai_wu_criteria(sigma, mat, opt) -> Float64

Tsai-Wu polynomial failure criterion.
Returns the margin of safety: sf - 1, where sf is the safety factor
at failure (positive = safe, negative = failed).
"""
function tsai_wu_criteria(sigma::Matrix{Float64}, mat::MaterialProperties, opt::String)
    size(sigma) == (3, 3) || throw(DimensionMismatch("stress tensor must be 3×3"))
    f1t, f1c, f2t, f2c, shear = _failure_strengths(mat, opt)
    σ1, σ2, τ = sigma[1, 1], sigma[2, 2], sigma[1, 2]
    f1, f11 = 1 / f1t - 1 / f1c, 1 / (f1t * f1c)
    f2, f22 = 1 / f2t - 1 / f2c, 1 / (f2t * f2c)
    f66 = 1 / shear^2
    f12 = -0.5 * sqrt(f11 * f22)
    linear = f1 * σ1 + f2 * σ2
    quadratic = f11 * σ1^2 + f22 * σ2^2 + f66 * τ^2 + 2.0f12 * σ1 * σ2
    iszero(linear) && iszero(quadratic) && return Inf
    roots = quadratic_roots(quadratic, linear, -1.0)
    valid_roots = filter(root -> isreal(root) && real(root) >= 0, roots)
    safety_factor = maximum(real.(valid_roots); init = 0.0)
    return safety_factor - 1
end

"""
    tsai_hill_criteria(sigma, mat) -> Float64

Tsai-Hill failure criterion. Returns 1 - FI (positive = safe, 0 = onset of failure).
"""
function tsai_hill_criteria(sigma::Matrix{Float64}, mat::MaterialProperties)
    size(sigma) == (3, 3) || throw(DimensionMismatch("stress tensor must be 3×3"))
    σ1, σ2, τ = sigma[1, 1], sigma[2, 2], sigma[1, 2]
    f1 = σ1 >= 0 ? mat.axial_trac : mat.axial_comp
    f2 = σ2 >= 0 ? mat.trans_trac : mat.trans_comp
    index = (σ1 / f1)^2 + (σ2 / f2)^2 - σ1 * σ2 / f1^2 + (τ / mat.shear)^2
    return 1 - index
end
