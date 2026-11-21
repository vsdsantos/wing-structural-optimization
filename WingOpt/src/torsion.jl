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
        throw(
            DimensionMismatch(
                "one thickness and segment length are required per wall segment",
            ),
        )
    all(>(0), t_vec) || throw(DomainError(t_vec, "wall thicknesses must be positive"))
    all(>=(0), s_vec) ||
        throw(DomainError(s_vec, "wall segment lengths must be nonnegative"))
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
    length(lim) == 2 || throw(DimensionMismatch("lim must contain the start and end span"))
    Am > 0 && G > 0 ||
        throw(DomainError((Am, G), "area and shear modulus must be positive"))
    wall_integral = integral_linha_media(thick_vec, s_vec)
    dphi_dx = _ -> torque * wall_integral / (4G * Am^2)
    return (dphi_dx(0.0) * (lim[2] - lim[1]), dphi_dx)
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
    length(torque_vec) == length(area_vec) == 2 ||
        throw(DimensionMismatch("two cell torques and areas are required"))
    length(t_cells) == length(s_cells) == 3 ||
        throw(DimensionMismatch("three wall groups are required"))
    all(>(0), area_vec) || throw(DomainError(area_vec, "cell areas must be positive"))
    integrals = [integral_linha_media(t_cells[i], s_cells[i]) for i in 1:3]
    matrix = [
        rel_q -rel_q 0.0
        2area_vec[1] 2area_vec[2] 0.0
        integrals[1] / area_vec[1] -integrals[2] / area_vec[2] 2integrals[3] *
            sum(1 ./ area_vec)
    ]
    return matrix \ [0.0, sum(torque_vec), 0.0]
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
    n = length(stations)
    n == length(props) == length(Mt) ||
        throw(
            DimensionMismatch(
                "stations, section properties, and torque must have equal lengths",
            ),
        )
    n >= 2 || throw(ArgumentError("at least two spanwise stations are required"))
    mat.G > 0 || throw(DomainError(mat.G, "shear modulus must be positive"))
    span = [station.span for station in stations]
    issorted(span) || throw(ArgumentError("stations must be ordered from root to tip"))
    areas = [prop.area_mid for prop in props]
    all(>(0), areas) ||
        throw(DomainError(areas, "mid-line enclosed areas must be positive"))
    all(>(0), getfield.(stations, :thickness)) ||
        throw(
            DomainError(
                getfield.(stations, :thickness),
                "wall thicknesses must be positive",
            ),
        )
    rate = Mt .* (getfield.(props, :height) .+ getfield.(props, :width)) ./
           (2 .* areas .^ 2 .* mat.G .* getfield.(stations, :thickness))
    return cumtrapz(span, rate)
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
    n >= 3 || throw(ArgumentError("at least three spanwise stations are required"))
    n == length(props) == length(load.R.M) == length(load.L.M) ||
        throw(
            DimensionMismatch(
                "stations, properties, and bending-moment arrays must have equal lengths",
            ),
        )
    span = [station.span for station in stations]
    increments = diff(span)
    all(>(0), increments) ||
        throw(ArgumentError("stations must be ordered from root to tip"))
    h = abs(increments[1])
    all(isapprox.(increments, h; rtol = 1e-8, atol = 1e-12)) ||
        throw(ArgumentError("wing_displacement requires uniformly spaced stations"))
    inertia = [prop.Iz0 for prop in props]
    all(>(0), inertia) || throw(DomainError(inertia, "bending inertias must be positive"))
    mat.axial_E > 0 || throw(DomainError(mat.axial_E, "Young's modulus must be positive"))

    matrix = zeros(n, n)
    matrix[1, 2] = 2
    matrix[2, 2:3] .= [-2, 1]
    for i in 3:(n - 1)
        matrix[i, (i - 1):(i + 1)] .= [1, -2, 1]
    end
    interior = matrix[2:(n - 1), 2:(n - 1)]
    function half_deflection(moment)
        curvature_rhs = -h^2 .* moment ./ (mat.axial_E .* inertia)
        slope = interior \ curvature_rhs[2:(n - 1)]
        return cumsum([0.0; slope; 0.0])
    end
    return (R = half_deflection(load.R.M), L = half_deflection(load.L.M))
end
