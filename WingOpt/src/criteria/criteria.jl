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
    1 <= idx <= length(load_half.s) ||
        throw(BoundsError(load_half.s, idx))
    length(load_half.V) == length(load_half.M) == length(load_half.Mt) ==
        length(load_half.F) ||
        throw(DimensionMismatch("all internal-load arrays must have equal lengths"))
    t_ref, t_box = station.reinforcer, station.thickness
    t_ref > 0 && t_box > 0 ||
        throw(
            DomainError(
                (t_ref, t_box),
                "wall and reinforcer thicknesses must be positive",
            ),
        )
    width, height = prop.width, prop.height
    width > 0 && height > 0 && prop.area > 0 && prop.Iz0 > 0 && prop.area_mid > 0 ||
        throw(DomainError(prop, "section dimensions, area, and inertia must be positive"))

    V, M, F, Mt = load_half.V[idx], load_half.M[idx], load_half.F[idx], load_half.Mt[idx]
    centroid_x = prop.centroid[1]
    x_web1 = station.limits[1] * station.chord - station.chord / 4
    x_web2 = station.limits[2] * station.chord - station.chord / 4
    y_vertices = prop.vertices[:, 2] .- prop.centroid[2]
    neutral_axis = iszero(M) ? 0.0 : F * prop.Iz0 / (M * prop.area)
    compression_points = y_vertices[y_vertices .>= neutral_axis]
    tension_points = y_vertices[y_vertices .< neutral_axis]
    y_comp = isempty(compression_points) ? maximum(y_vertices) : maximum(compression_points)
    y_tens = isempty(tension_points) ? minimum(y_vertices) : minimum(tension_points)

    boom_area = t_ref^2
    boom_centroid = height / 2
    boom_inertia = 2boom_area * ((height - boom_centroid)^2 + boom_centroid^2)
    sigma_x = iszero(boom_inertia) ? zeros(4) :
              -M .* [y_comp, y_comp, y_tens, y_tens] ./ boom_inertia

    denominator = iszero(y_tens) ? copysign(eps(Float64), y_tens) : y_tens
    ratio_ct = y_comp / denominator
    ratio_tc = iszero(y_comp) ? 0.0 : y_tens / y_comp
    boom_shear_area = t_box * (height * (2 + ratio_ct) / 6 + width / 2)
    boom_tension_area = t_box * (height * (2 + ratio_tc) / 6 + width / 2)
    boom_shear_area > 0 && boom_tension_area > 0 ||
        throw(
            DomainError(
                (boom_shear_area, boom_tension_area),
                "shear-flow idealization produced nonpositive boom areas",
            ),
        )
    shear_centroid = boom_shear_area * height / (boom_shear_area + boom_tension_area)
    shear_inertia = 2boom_shear_area * (height - shear_centroid)^2 +
                    2boom_tension_area * shear_centroid^2
    shear_inertia > 0 ||
        throw(DomainError(shear_inertia, "shear-flow inertia must be positive"))

    qb12 = 0.0
    qb23 = -(V / shear_inertia) * (boom_shear_area * y_comp)
    qb34 = qb23 - (V / shear_inertia) * (boom_tension_area * y_tens)
    qb41 = qb34
    area23 = (centroid_x - x_web1) * height / 4
    area34 = height * width / 4
    area41 = (x_web2 - centroid_x) * height / 4
    q0 = -2 * (area23 * qb23 + area34 * qb34 + area41 * qb41) / (2 * width * height)
    q_t = Mt / (2 * width * height)
    shear_flows = [q0 + q + q_t for q in (qb12, qb23, qb34, qb41)]
    shear_stresses = shear_flows ./ t_box

    box_material = material_properties("balsa")
    reinforcer_material = material_properties("freijo")
    margins = Float64[]
    for (boom_idx, material) in ((2, reinforcer_material), (3, reinforcer_material),
                                (2, box_material), (4, box_material),
                                (1, box_material), (3, box_material))
        sigma = zeros(3, 3)
        if length(margins) <= 1
            sigma[1, 1] = sigma_x[boom_idx]
            fiber = "axial"
        else
            flow_idx = length(margins) - 1
            sigma[1, 2] = sigma[2, 1] = shear_stresses[flow_idx]
            fiber = length(margins) >= 4 ? "trans" : "axial"
        end
        push!(margins, max_tension_criteria(sigma, material, fiber))
    end
    # Apply the Tsai-Wu check to the same two boom and four wall locations.
    for (location, (boom_idx, material)) in enumerate((
        (2, reinforcer_material), (3, reinforcer_material),
        (2, box_material), (4, box_material), (1, box_material), (3, box_material),
    ))
        sigma = zeros(3, 3)
        if location <= 2
            sigma[1, 1] = sigma_x[boom_idx]
            fiber = "axial"
        else
            flow_idx = location - 2
            sigma[1, 2] = sigma[2, 1] = shear_stresses[flow_idx]
            fiber = location >= 5 ? "trans" : "axial"
        end
        push!(margins, tsai_wu_criteria(sigma, material, fiber))
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
    n = length(stations)
    n == length(props) ||
        throw(DimensionMismatch("one section property record is required per station"))
    n >= 2 || throw(ArgumentError("at least two stations are required to integrate mass"))
    span = [station.span for station in stations]
    increments = diff(span)
    all(>(0), increments) ||
        throw(ArgumentError("stations must be ordered from root to tip"))

    box_area = [
        2 * prop.width * station.thickness + 2 * prop.height * station.thickness
        for (station, prop) in zip(stations, props)
    ]
    reinforcer_area = [2 * station.reinforcer^2 for station in stations]
    all(>=(0), box_area) && all(>=(0), reinforcer_area) ||
        throw(
            DomainError(
                (box_area, reinforcer_area),
                "cross-sectional areas must be nonnegative",
            ),
        )
    dm = [
        increments[i] * (
            (box_area[i] + box_area[i + 1]) * mat_box.dens +
            (reinforcer_area[i] + reinforcer_area[i + 1]) * mat_reinforcer.dens
        ) / 2
        for i in 1:(n - 1)
    ]
    return (sum(dm), dm)
end

"""
    report_45(stations, props, wing, mat; design_speed) -> Float64

Compute the wing flexibility factor F from ANC-45 (Report 45).

    F = Σ (θᵢ/M) · Cᵢ² · Δsᵢ

where θᵢ is the torsion angle per unit torque at station i (from the aileron
cutout span outward), Cᵢ is the local chord in feet, and Δsᵢ is the
span increment in feet.

`design_speed` is the design dive speed in m/s. The criterion is F ≤ 200/Vd²,
where Vd is converted to mph. Returns F / (200/Vd²): values < 1 pass.
"""
function report_45(
    stations::Vector{Station},
    props::Vector{SectionProperties},
    wing::WingGeometry,
    mat::MaterialProperties,
    ;
    design_speed::Float64,
)
    design_speed > 0 || throw(DomainError(design_speed, "design speed must be positive"))
    length(stations) == length(props) ||
        throw(DimensionMismatch("one section property record is required per station"))
    length(wing.span_transitions) >= 2 ||
        throw(ArgumentError("wing must include a start and end span transition"))
    span = [station.span for station in stations]
    issorted(span) || throw(ArgumentError("stations must be ordered from root to tip"))
    span[1] <= wing.span_transitions[2] <= span[end] ||
        throw(
            DomainError(
                wing.span_transitions[2],
                "aileron cutout must lie within the stations",
            ),
        )
    start_index = argmin(abs.(span .- wing.span_transitions[2]))
    twist = wing_torsion(stations, props, ones(length(stations)), mat)
    remaining_span = span[start_index:end]
    tip_span = wing.span_transitions[end]
    tip_span >= remaining_span[end] ||
        throw(DomainError(tip_span, "wing tip cannot be inboard of the analyzed stations"))
    trailing_width = tip_span - remaining_span[end]
    widths = [diff(remaining_span); trailing_width]
    edges = [remaining_span .- widths ./ 2; tip_span]
    delta_span = diff(edges) .* 3.28084
    chord = [station.chord for station in stations[start_index:end]] .* 3.28084
    design_speed_mph = design_speed * 2.23694
    flexibility = sum(
        (twist[start_index:end] ./ 0.73756) .* chord .^ 2 .* delta_span
    )
    return flexibility / (200 / design_speed_mph^2)
end
