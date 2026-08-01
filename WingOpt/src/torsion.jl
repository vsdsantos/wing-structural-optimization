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

Compute `∮ ds/t` around a closed thin-walled section, given vectors of
wall thickness `t_vec` and segment lengths `s_vec` (both length N).
"""
function integral_linha_media(t_vec::AbstractVector, s_vec::AbstractVector)
    error("Not implemented")
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
    torque::Float64,
    thick_vec::AbstractVector,
    s_vec::AbstractVector,
    Am::Float64,
    G::Float64,
    lim::AbstractVector,
)
    error("Not implemented")
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
    rel_q::Float64,
    t_cells::AbstractVector,
    s_cells::AbstractVector,
)
    error("Not implemented")
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
    error("Not implemented")
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
