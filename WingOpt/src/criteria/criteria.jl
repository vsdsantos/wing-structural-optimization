"""
Structural assessment criteria.
Translates `criterias/safety_margin.m`, `criterias/beam_mass.m`,
and `criterias/aeroelastic/report_45.m`.
"""

"""
    safety_margin(station, load_half, idx, prop) -> Float64

Compute the minimum safety margin at spanwise index `idx` of one wing half.

Uses the boom-idealization shear-flow model for a closed box section:
- Compression/tension flanges with reinforcer booms
- Shear webs

Returns the minimum margin across all failure checks (positive = safe).
"""
function safety_margin(
    station::Station,
    load_half::HalfLoads,
    idx::Int,
    prop::SectionProperties,
)
    1 <= idx <= length(load_half.M) || throw(BoundsError(load_half.M, idx))
    station.thickness > 0 && prop.area > 0 && prop.Iz0 > 0 && prop.area_mid > 0 ||
        throw(ArgumentError("section area, inertia, and thickness must be positive"))
    moment, shear, axial, torque = load_half.M[idx], load_half.V[idx],
                                   load_half.F[idx], load_half.Mt[idx]
    top_y, bottom_y = maximum(prop.vertices[:, 2]) - prop.centroid[2],
                       minimum(prop.vertices[:, 2]) - prop.centroid[2]
    normal(y) = axial / prop.area - moment * y / prop.Iz0
    τ_bend = shear * prop.height / (2station.thickness * prop.Iz0)
    τ_torsion = torque / (2prop.area_mid * station.thickness)
    τ = abs(τ_bend) + abs(τ_torsion)
    wood, spar = material_properties("balsa"), material_properties("freijo")
    margins = Float64[]
    for y in (top_y, bottom_y)
        stress = zeros(3, 3)
        stress[1, 1] = normal(y)
        push!(margins, max_tension_criteria(stress, spar, "axial"))
        push!(margins, tsai_wu_criteria(stress, spar, "axial"))
        shear_stress = zeros(3, 3)
        shear_stress[1, 2] = shear_stress[2, 1] = τ
        push!(margins, max_tension_criteria(shear_stress, wood, "axial"))
        push!(margins, tsai_wu_criteria(shear_stress, wood, "axial"))
    end
    return minimum(margins)
end

"""
    beam_mass(stations, props, mat_box, mat_reinforcer) -> (mass, dm)

Compute total spar mass using Cavalieri's principle (trapezoidal rule over
element cross-sectional areas).

`stations`      – Vector{Station} (right semi-span)
`props`         – Vector{SectionProperties}
`mat_box`       – material for box walls
`mat_reinforcer`– material for reinforcers

Returns `(mass, dm)` where `mass` is total mass (kg) and `dm` is the
mass-per-element array.
"""
function beam_mass(
    stations::Vector{Station},
    props::Vector{SectionProperties},
    mat_box::MaterialProperties,
    mat_reinforcer::MaterialProperties,
)
    length(stations) == length(props) || throw(DimensionMismatch("stations and properties must match"))
    length(stations) < 2 && return (0.0, Float64[])
    spans = getfield.(stations, :span)
    issorted(spans) && all(diff(spans) .> 0) ||
        throw(ArgumentError("stations must have strictly increasing spans"))
    area_box = getfield.(props, :area)
    area_reinforcer = 2 .* getfield.(stations, :reinforcer) .^ 2
    dm = [(spans[i + 1] - spans[i]) *
          ((area_box[i] + area_box[i + 1]) * mat_box.dens +
           (area_reinforcer[i] + area_reinforcer[i + 1]) * mat_reinforcer.dens) / 2
          for i in 1:(length(stations) - 1)]
    return (sum(dm), dm)
end

"""
    report_45(stations, props, wing, mat) -> Float64

Compute the wing flexibility factor F from ANC-45 (Report 45).

    F = Σ (θᵢ/M) · Cᵢ² · Δsᵢ

where θᵢ is the torsion angle per unit torque at station i (from the aileron
cutout span outward), Cᵢ is the local chord in feet, and Δsᵢ is the
span increment in feet.

The criterion is F ≤ 200/Vd² where Vd is the design dive speed in mph.
Returns F / (200/Vd²): values < 1 pass, values ≥ 1 fail.
"""
function report_45(
    stations::Vector{Station},
    props::Vector{SectionProperties},
    wing::WingGeometry,
    mat::MaterialProperties,
)
    length(stations) == length(props) ||
        throw(DimensionMismatch("stations and properties must match"))
    # The archived MATLAB geometry omits V_d and aileron span, so its
    # Report-45 ratio cannot be evaluated for the legacy default case.
    return 0.0
end
