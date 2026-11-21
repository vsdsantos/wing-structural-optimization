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
    length(t_vec) == length(s_vec) || throw(DimensionMismatch("thickness and segment lengths must match"))
    all(>(0), t_vec) || throw(ArgumentError("wall thicknesses must be positive"))
    all(>=(0), s_vec) || throw(ArgumentError("segment lengths must be nonnegative"))
    return sum(s_vec ./ t_vec)
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
    Am > 0 && G > 0 || throw(ArgumentError("area and shear modulus must be positive"))
    length(lim) == 2 && lim[2] >= lim[1] || throw(ArgumentError("lim must be [start, end]"))
    compliance = integral_linha_media(thick_vec, s_vec) / (4G * Am^2)
    torque_at(x) = torque isa Function ? torque(x) : torque
    derivative = x -> torque_at(x) .* compliance
    samples = range(Float64(lim[1]), Float64(lim[2]); length = 201)
    total = trapz(collect(samples), derivative.(samples))
    return total, derivative
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
    length(torque_vec) == 2 && length(area_vec) == 2 && length(t_cells) == 3 && length(s_cells) == 3 ||
        throw(DimensionMismatch("bi-cellular torsion expects 2 torques, 2 areas, and 3 wall groups"))
    all(>(0), area_vec) || throw(ArgumentError("cell areas must be positive"))
    integrals = [integral_linha_media(t_cells[i], s_cells[i]) for i in 1:3]
    matrix = [rel_q 1.0 1.0;
              2area_vec[1] 2area_vec[2] 0.0;
              integrals[1] / area_vec[1] -integrals[2] / area_vec[2]
              2integrals[3] * (1 / area_vec[1] + 1 / area_vec[2])]
    rhs = [0.0, sum(torque_vec), 0.0]
    return matrix \ rhs
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
    length(stations) == length(props) == length(Mt) || throw(DimensionMismatch("stations, properties, and torque arrays must match"))
    isempty(stations) && return Float64[]
    span = getfield.(stations, :span)
    area = getfield.(props, :area_mid)
    all(>(0), area) || throw(ArgumentError("mid-line enclosed areas must be positive"))
    gradient = [Mt[i] * 2(props[i].width + props[i].height) /
                (4area[i]^2 * mat.G * stations[i].thickness) for i in eachindex(stations)]
    return cumtrapz(span, gradient)
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
    n = length(stations)
    n == length(props) || throw(DimensionMismatch("stations and properties must match"))
    n == length(load.R.M) && n == length(load.L.M) ||
        throw(DimensionMismatch("load arrays must match the station count"))
    n < 3 && return (R = zeros(n), L = zeros(n))
    h = stations[2].span - stations[1].span
    h > 0 || throw(ArgumentError("stations must have strictly increasing span coordinates"))
    function deflection(moment)
        matrix = zeros(n - 2, n - 2)
        rhs = zeros(n - 2)
        for i in 1:(n - 2)
            matrix[i, i] = -2
            i > 1 && (matrix[i, i - 1] = 1)
            i < n - 2 && (matrix[i, i + 1] = 1)
            rhs[i] = -h^2 * moment[i + 1] / (mat.axial_E * props[i + 1].Iz0)
        end
        interior = matrix \ rhs
        return cumsum(vcat(0.0, interior, 0.0))
    end
    return (R = deflection(load.R.M), L = deflection(load.L.M))
end
