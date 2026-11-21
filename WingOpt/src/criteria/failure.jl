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
function max_tension_criteria(
    sigma::AbstractMatrix{<:Real},
    mat::MaterialProperties,
    opt::String,
)
    size(sigma, 1) >= 2 && size(sigma, 2) >= 2 ||
        throw(DimensionMismatch("stress tensor must contain its in-plane 2×2 components"))
    strengths = _failure_strengths(mat, opt)
    all(>(0), strengths) && mat.shear > 0 ||
        throw(DomainError(strengths, "material strengths must be positive"))
    σ1, σ2, τ12 = sigma[1, 1], sigma[2, 2], sigma[1, 2]
    margins = Float64[]
    for (stress, tensile, compressive) in (
        (σ1, strengths[1], strengths[2]),
        (σ2, strengths[3], strengths[4]),
    )
        if stress > 0
            push!(margins, tensile / stress - 1)
        elseif stress < 0
            push!(margins, -compressive / stress - 1)
        end
    end
    !iszero(τ12) && push!(margins, mat.shear / abs(τ12) - 1)
    return isempty(margins) ? Inf : minimum(margins)
end

"""
    tsai_wu_criteria(sigma, mat, opt) -> Float64

Tsai-Wu polynomial failure criterion.
Returns the margin of safety: sf - 1, where sf is the safety factor
at failure (positive = safe, negative = failed).
"""
function tsai_wu_criteria(sigma::Matrix{Float64}, mat::MaterialProperties, opt::String)
    size(sigma, 1) >= 2 && size(sigma, 2) >= 2 ||
        throw(DimensionMismatch("stress tensor must contain its in-plane 2×2 components"))
    F1t, F1c, F2t, F2c = _failure_strengths(mat, opt)
    F6 = mat.shear
    all(>(0), (F1t, F1c, F2t, F2c, F6)) ||
        throw(DomainError((F1t, F1c, F2t, F2c, F6), "strengths must be positive"))

    σ1, σ2, τ12 = sigma[1, 1], sigma[2, 2], sigma[1, 2]
    f1 = 1 / F1t - 1 / F1c
    f11 = 1 / (F1t * F1c)
    f2 = 1 / F2t - 1 / F2c
    f22 = 1 / (F2t * F2c)
    f66 = 1 / F6^2
    f12 = -0.5 * sqrt(f11 * f22)
    linear = f1 * σ1 + f2 * σ2
    quadratic = f11 * σ1^2 + f22 * σ2^2 + f66 * τ12^2 + f12 * σ1 * σ2
    if iszero(linear) && iszero(quadratic)
        return Inf
    end
    r1, r2 = quadratic_roots(quadratic, linear, -1.0)
    factor = maximum((r1, r2))
    return factor - 1
end

"""
    tsai_hill_criteria(sigma, mat) -> Float64

Tsai-Hill failure criterion. Returns the Tsai-Hill failure index (values below 1 are safe).
"""
function tsai_hill_criteria(sigma::Matrix{Float64}, mat::MaterialProperties)
    size(sigma, 1) >= 2 && size(sigma, 2) >= 2 ||
        throw(DimensionMismatch("stress tensor must contain its in-plane 2×2 components"))
    σ1, σ2, τ12 = sigma[1, 1], sigma[2, 2], sigma[1, 2]
    F1 = σ1 >= 0 ? mat.axial_trac : mat.axial_comp
    F2 = σ2 >= 0 ? mat.trans_trac : mat.trans_comp
    all(>(0), (F1, F2, mat.shear)) ||
        throw(DomainError((F1, F2, mat.shear), "strengths must be positive"))
    return (σ1 / F1)^2 + (σ2 / F2)^2 - σ1 * σ2 / F1^2 + (τ12 / mat.shear)^2
end

function _failure_strengths(mat::MaterialProperties, opt::String)
    if opt == "axial"
        return (mat.axial_trac, mat.axial_comp, mat.trans_trac, mat.trans_comp)
    elseif opt == "trans"
        return (mat.trans_trac, mat.trans_comp, mat.axial_trac, mat.axial_comp)
    end
    throw(ArgumentError("opt must be \"axial\" or \"trans\""))
end
