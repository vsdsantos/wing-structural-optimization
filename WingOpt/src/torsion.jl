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
    length(t_vec) == length(s_vec) ||
        error("t_vec and s_vec must have equal length")
    return sum(s_vec[i] / t_vec[i] for i in eachindex(t_vec))
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
function ang_torcao_parede_fina(torque::Float64, thick_vec::AbstractVector,
                                  s_vec::AbstractVector, Am::Float64,
                                  G::Float64, lim::AbstractVector)
    int_val  = integral_linha_media(thick_vec, s_vec)
    dphi_dx  = x -> torque * int_val / (4 * G * Am^2)
    phi_total = dphi_dx(0.0) * (lim[2] - lim[1])
    return phi_total, dphi_dx
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
function bicelular(torque_vec::AbstractVector, area_vec::AbstractVector,
                   rel_q::Float64,
                   t_cells::AbstractVector, s_cells::AbstractVector)
    length(t_cells) == 3 && length(s_cells) == 3 ||
        error("t_cells and s_cells must each have length 3")
    int = [integral_linha_media(t_cells[i], s_cells[i]) for i in 1:3]
    A = [rel_q                              0.0                          0.0;
         2*area_vec[1]                      2*area_vec[2]                0.0;
         int[1]/area_vec[1]  -int[2]/area_vec[2]  2*int[3]*sum(1 ./ area_vec)]
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
function wing_torsion(stations::Vector{Station},
                      props::Vector{SectionProperties},
                      Mt::AbstractVector, mat::MaterialProperties)
    sp = [s.span for s in stations]
    Am = [p.area_mid for p in props]
    tc = [s.thickness for s in stations]
    h  = [p.height   for p in props]
    w  = [p.width    for p in props]

    dphidz = Mt .* 2 .* (h .+ w) ./ (4 .* Am.^2 .* mat.G .* tc)
    return cumtrapz(sp, dphidz)
end

# ---------------------------------------------------------------------------
# Wing-level displacement (from criterias/displacement/wing_displacement.m)
# ---------------------------------------------------------------------------

"""
    wing_displacement(stations, props, load, mat) -> NamedTuple

Finite-difference beam bending displacement for both wing halves.

Returns `(R = right_deflection, L = left_deflection)` as vectors over `stations`.
"""
function wing_displacement(stations::Vector{Station},
                            props::Vector{SectionProperties},
                            load::InternalLoads, mat::MaterialProperties)
    s = [st.span for st in stations]
    n = length(s)
    h = abs(s[2] - s[1])
    E = mat.axial_E
    I_vals = [p.Iz0 for p in props]

    K_R = -h^2 .* load.R.M ./ (E .* I_vals)
    K_L = -h^2 .* load.L.M ./ (E .* I_vals)

    # Build tri-diagonal stiffness matrix (clamped at root, free at tip)
    A = zeros(n, n)
    A[1, 2] = 2.0
    n >= 2 && (A[2, 2:3] = [-2.0, 1.0])
    for i in 3:n-1
        A[i, i-1:i+1] = [1.0, -2.0, 1.0]
    end
    # Remove boundary rows/cols (pinned at root and tip)
    idx = 2:n-1
    A_inner = A[idx, idx]
    K_R_inner = K_R[idx]
    K_L_inner = K_L[idx]

    w_r = A_inner \ K_R_inner
    w_l = A_inner \ K_L_inner

    right = vcat([0.0], w_r, [0.0])
    left  = vcat([0.0], w_l, [0.0])
    return (R = cumsum(right), L = cumsum(left))
end
