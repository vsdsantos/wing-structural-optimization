"""
Torsion calculations for thin-walled beam sections.
Translates `torsion/`: integral_linha_media.m, ang_torcao_parede_fina.m,
bicelular.m; and `criterias/displacement/wing_torsion.m`.
"""

# ---------------------------------------------------------------------------
# Closed thin-wall section (Bredt's formula)
# ---------------------------------------------------------------------------

"""
    integral_linha_media(t_vec, s_vec) -> Float64

Compute the MATLAB implementation's `sum(tᵢ * sᵢ)` around a closed
thin-walled section. This preserves its current input and numerical convention.
"""
function integral_linha_media(t_vec::AbstractVector, s_vec::AbstractVector)
    length(t_vec) == length(s_vec) ||
        throw(DimensionMismatch("thickness and segment-length vectors must match"))
    isempty(t_vec) && return 0.0
    all(isfinite, t_vec) && all(isfinite, s_vec) ||
        throw(ArgumentError("thicknesses and segment lengths must be finite"))
    all(>(0), t_vec) && all(>=(0), s_vec) ||
        throw(ArgumentError("thicknesses must be positive and segment lengths nonnegative"))
    return sum(Float64(t) * Float64(s) for (t, s) in zip(t_vec, s_vec))
end

"""
    ang_torcao_parede_fina(torque, thick_vec, s_vec, Am, G, lim) -> (phi_total, dphi_dx)

Torsion angle for a single-cell thin-walled section (Bredt–Batho).

- `torque`    – applied torque (N·m); constant scalar
- `thick_vec` – vector of wall thicknesses for each segment
- `s_vec`     – vector of segment arc-lengths
- `Am`        – enclosed mid-line area (m²)
- `G`         – shear modulus (Pa)
- `lim`       – integration limits [x_start, x_end] along the span

Returns `(phi_total, dphi_dx)` where `dphi_dx` is a function `x → dφ/dx`.
"""
function ang_torcao_parede_fina(
    torque::Union{Real,Function},
    thick_vec::AbstractVector,
    s_vec::AbstractVector,
    Am::Real,
    G::Real,
    lim::AbstractVector,
)
    length(lim) == 2 || throw(ArgumentError("lim must contain two span limits"))
    isfinite(Am) && Am > 0 || throw(ArgumentError("Am must be positive and finite"))
    isfinite(G) && G > 0 || throw(ArgumentError("G must be positive and finite"))
    t_integral = integral_linha_media(thick_vec, s_vec)
    torque_fn = torque isa Function ? torque : x -> zero.(x) .+ torque
    dphi_dx = x -> torque_fn(x) .* (t_integral / (4 * G * Am^2))
    phi_total = _simpson_integral(dphi_dx, Float64(lim[1]), Float64(lim[2]))
    return phi_total, dphi_dx
end

function _simpson_integral(f, a::Float64, b::Float64)
    isfinite(a) && isfinite(b) || throw(ArgumentError("integration limits must be finite"))
    a == b && return 0.0
    n = 256
    h = (b - a) / n
    total = Float64(f(a)) + Float64(f(b))
    for i = 1:(n-1)
        total += (isodd(i) ? 4 : 2) * Float64(f(a + i * h))
    end
    return total * h / 3
end

# ---------------------------------------------------------------------------
# Bi-cellular section
# ---------------------------------------------------------------------------

"""
    bicelular(torque_vec, area_vec, rel_q, t_cells, s_cells) -> Vector{Float64}

Compute shear flows in a bi-cellular thin-walled section.

- `torque_vec` – [Mt1, Mt2] torques in each cell (N·m)
- `area_vec`   – [A1, A2] enclosed areas (m²)
- `rel_q`      – shear-flow compatibility parameter
- `t_cells`    – length-3 vector of thickness vectors (one per wall group)
- `s_cells`    – length-3 vector of segment-length vectors

Returns a 3-element vector of shear flows.
"""
function bicelular(
    torque_vec::AbstractVector,
    area_vec::AbstractVector,
    rel_q::AbstractVector,
    t_cells::AbstractVector,
    s_cells::AbstractVector,
)
    length(torque_vec) == 2 && length(area_vec) == 2 ||
        throw(DimensionMismatch("the bi-cell model requires two torques and areas"))
    length(rel_q) == length(t_cells) == length(s_cells) == 3 ||
        throw(DimensionMismatch("rel_q, t_cells and s_cells must each have length 3"))
    all(isfinite, area_vec) && all(>(0), area_vec) ||
        throw(ArgumentError("cell areas must be positive and finite"))
    all(isfinite, torque_vec) && all(isfinite, rel_q) ||
        throw(ArgumentError("torques and compatibility coefficients must be finite"))
    line_integrals =
        [integral_linha_media(t_cells[i], s_cells[i]) for i in eachindex(t_cells)]
    A = [
        Float64.(rel_q)'
        2 * area_vec[1] 2 * area_vec[2] 0
        line_integrals[1] / area_vec[1] -line_integrals[2] / area_vec[2] 2 * line_integrals[3] * sum(1 ./ area_vec)
    ]
    B = [0.0, sum(torque_vec), 0.0]
    return A \ B
end

# ---------------------------------------------------------------------------
# Wing-level torsion (from criterias/displacement/wing_torsion.m)
# ---------------------------------------------------------------------------

"""
    wing_torsion(stations, props, Mt, mat) -> Vector{Float64}

Cumulative torsion angle along the span.

- `stations` – Vector{Station}
- `props`    – Vector{SectionProperties}
- `Mt`       – torsional moment at each station (N·m)
- `mat`      – MaterialProperties (uses `.G`)
"""
function wing_torsion(
    stations::Vector{Station},
    props::Vector{SectionProperties},
    Mt::AbstractVector,
    mat::MaterialProperties,
)
    length(stations) == length(props) == length(Mt) || throw(
        DimensionMismatch("stations, properties and torque must have matching lengths"),
    )
    isempty(stations) && return Float64[]
    all(isfinite, Mt) || throw(ArgumentError("torques must be finite"))
    spans = getfield.(stations, :span)
    length(spans) < 2 ||
        all(diff(spans) .> 0) ||
        throw(ArgumentError("station spans must be strictly increasing"))
    G = mat.G
    isfinite(G) && G > 0 ||
        throw(ArgumentError("shear modulus must be positive and finite"))
    dphi_dx = similar(Float64.(Mt))
    for i in eachindex(stations)
        Am = props[i].area_mid
        t = stations[i].thickness
        h = props[i].height
        w = props[i].width
        isfinite(Am) && Am > 0 ||
            throw(ArgumentError("mid-line area must be positive and finite"))
        isfinite(t) && t > 0 ||
            throw(ArgumentError("wall thickness must be positive and finite"))
        dphi_dx[i] = Mt[i] * 2 * (h + w) / (4 * Am^2 * G * t)
    end
    return cumtrapz(spans, dphi_dx)
end

# ---------------------------------------------------------------------------
# Wing-level displacement (from criterias/displacement/wing_displacement.m)
# ---------------------------------------------------------------------------

"""
    wing_displacement(stations, props, load, mat) -> NamedTuple

Finite-difference beam bending displacement for both wing halves.

Returns `(R = right_deflection, L = left_deflection)` as vectors over `stations`.
"""
function wing_displacement(
    stations::Vector{Station},
    props::Vector{SectionProperties},
    load::InternalLoads,
    mat::MaterialProperties,
)
    error("Not implemented")
end
