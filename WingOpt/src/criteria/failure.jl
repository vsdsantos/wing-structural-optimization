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
function max_tension_criteria(sigma::Matrix{Float64},
                               mat::MaterialProperties, opt::String)
    error("Not implemented")
end

"""
    tsai_wu_criteria(sigma, mat, opt) -> Float64

Tsai-Wu polynomial failure criterion.
Returns the margin of safety: sf - 1, where sf is the safety factor
at failure (positive = safe, negative = failed).
"""
function tsai_wu_criteria(sigma::Matrix{Float64},
                           mat::MaterialProperties, opt::String)
    error("Not implemented")
end

"""
    tsai_hill_criteria(sigma, mat) -> Float64

Tsai-Hill failure criterion. Returns 1 - FI (positive = safe, 0 = onset of failure).
"""
function tsai_hill_criteria(sigma::Matrix{Float64}, mat::MaterialProperties)
    error("Not implemented")
end
