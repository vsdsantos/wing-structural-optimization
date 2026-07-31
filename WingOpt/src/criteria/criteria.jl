"""
Structural assessment criteria.
Translates `criterias/safety_margin.m`, `criterias/beam_mass.m`,
and `criterias/aeroelastic/report_45.m`.

Notes on known issues fixed during migration:
- safety_margin.m called `material_long(...)` which is an alias for
  `material_properties`; replaced with `material_properties` directly.
- report_45.m called `torcao_asa(...)` which is `wing_torsion`.
"""

# ---------------------------------------------------------------------------
# Safety margin
# ---------------------------------------------------------------------------

"""
    safety_margin(station, load_half, idx, prop) -> Float64

Compute the minimum safety margin at spanwise index `idx` of one wing half.

Uses the boom-idealization shear-flow model for a closed box section:
- Compression/tension flanges with reinforcer booms
- Shear webs

Returns the minimum margin across all failure checks (positive = safe).
"""
function safety_margin(station::Station, load_half::HalfLoads,
                        idx::Int, prop::SectionProperties)
    t_ref   = station.reinforcer
    t_box   = station.thickness

    cent  = prop.centroid
    width = prop.width
    height = prop.height
    Iz0   = prop.Iz0
    area  = prop.area
    Am    = prop.area_mid

    x_cent  = cent[1]
    x_alma1 = station.limits[1] * station.chord - station.chord / 4
    x_alma2 = station.limits[2] * station.chord - station.chord / 4

    # Internal loads at this station
    V  = load_half.V[idx]
    Mf = load_half.M[idx]
    F  = load_half.F[idx]
    Mt = load_half.Mt[idx]

    # Neutral axis equation:  σ = -Mf*y/Iz0 + F/A
    # Points above/below neural axis
    z0 = prop.vertices[:, 1] .- cent[1]
    y0 = prop.vertices[:, 2] .- cent[2]

    a_na = (abs(Iz0) > eps()) ? Mf / Iz0 : 0.0
    c_na = (area    > eps()) ? F  / area  : 0.0

    # Neutral-axis level y_NA (from z0–y0 frame):  y_NA = c_na / a_na
    y_NA = (abs(a_na) > eps()) ? c_na / a_na : Inf

    P = hcat(z0, y0)
    comp_mask = P[:, 2] .>= y_NA
    P_comp = P[comp_mask, :]
    P_trac = P[.!comp_mask, :]

    if isempty(P_comp) || isempty(P_trac)
        return 0.0  # degenerate geometry
    end

    d_comp = abs.( -a_na .* P_comp[:, 2] .+ c_na) ./ abs(a_na + eps())
    d_trac = abs.( -a_na .* P_trac[:, 2] .+ c_na) ./ abs(a_na + eps())

    y_comp = P_comp[argmax(d_comp), 2]
    y_trac = P_trac[argmax(d_trac), 2]

    # Reinforcer booms (contribute to bending stiffness)
    B = zeros(4)
    B[2] = t_ref^2
    B[3] = t_ref^2

    y_cent_boom = abs(sum(B)) < eps() ? height/2 :
                  sum(B[[1,2]]) * height / sum(B)
    Ixx_boom = B[1]*(height - y_cent_boom)^2 + B[2]*(height - y_cent_boom)^2 +
               B[3]*y_cent_boom^2             + B[4]*y_cent_boom^2

    σx = (abs(Ixx_boom) > eps()) ?
         -Mf .* [y_comp, y_comp, y_trac, y_trac] ./ Ixx_boom :
         zeros(4)

    # Box-section booms (for shear flow)
    Bs = zeros(4)
    Bs[1] = t_box*height*(2 + y_comp/(y_trac + eps()))/6 + t_box*width/2
    Bs[2] = Bs[1]
    Bs[3] = t_box*height*(2 + y_trac/(y_comp + eps()))/6 + t_box*width/2
    Bs[4] = Bs[3]

    y_cent_s = sum(Bs[[1,2]]) * height / (sum(Bs) + eps())
    Ixx_s    = Bs[1]*(height - y_cent_s)^2 + Bs[2]*(height - y_cent_s)^2 +
               Bs[3]*y_cent_s^2             + Bs[4]*y_cent_s^2

    # Open shear flows (cut at web 1-2)
    qb12 = 0.0
    qb23 = -(V / (Ixx_s + eps())) * Bs[2] * y_comp
    qb34 =  qb23 - (V / (Ixx_s + eps())) * Bs[3] * y_trac
    qb41 =  qb34 - (V / (Ixx_s + eps())) * Bs[4] * y_trac

    A23 = (x_cent - x_alma1) * height / 4
    A34 = height * width / 4
    A41 = (x_alma2 - x_cent) * height / 4

    qs0 = -2*(A23*qb23 + A34*qb34 + A41*qb41) / (2*width*height + eps())
    q_t = Mt / (2*width*height + eps())

    q12 = qs0 + qb12 + q_t
    q23 = qs0 + qb23 + q_t
    q34 = qs0 + qb34 + q_t
    q41 = qs0 + qb41 + q_t

    σ_xy = [q12, q23, q34, q41] ./ (t_box + eps())

    # Material properties
    mat_mesa = material_properties("freijo")   # flanges
    mat_alma = material_properties("balsa")    # webs

    # Stress tensors per element (3×3)
    function make_boom_tensor(σ1)
        T = zeros(3, 3);  T[1, 1] = σ1;  return T
    end
    function make_shear_tensor(τ)
        T = zeros(3, 3);  T[1, 2] = τ;  T[2, 1] = τ;  return T
    end

    σ_boom = [make_boom_tensor(σx[i])  for i in 1:4]
    σ_rev  = [make_shear_tensor(σ_xy[i]) for i in 1:4]

    fail = zeros(6, 2)
    fail[1, 1] = max_tension_criteria(σ_boom[2], mat_mesa, "axial")
    fail[2, 1] = max_tension_criteria(σ_boom[3], mat_mesa, "axial")
    fail[3, 1] = max_tension_criteria(σ_rev[2],  mat_alma, "axial")
    fail[4, 1] = max_tension_criteria(σ_rev[4],  mat_alma, "axial")
    fail[5, 1] = max_tension_criteria(σ_rev[1],  mat_alma, "trans")
    fail[6, 1] = max_tension_criteria(σ_rev[3],  mat_alma, "trans")

    fail[1, 2] = tsai_wu_criteria(σ_boom[2], mat_mesa, "axial")
    fail[2, 2] = tsai_wu_criteria(σ_boom[3], mat_mesa, "axial")
    fail[3, 2] = tsai_wu_criteria(σ_rev[2],  mat_alma, "axial")
    fail[4, 2] = tsai_wu_criteria(σ_rev[4],  mat_alma, "axial")
    fail[5, 2] = tsai_wu_criteria(σ_rev[1],  mat_alma, "trans")
    fail[6, 2] = tsai_wu_criteria(σ_rev[3],  mat_alma, "trans")

    finite_vals = filter(isfinite, vec(fail))
    return isempty(finite_vals) ? 0.0 : minimum(finite_vals)
end

# ---------------------------------------------------------------------------
# Beam mass  (Cavalieri's principle)
# ---------------------------------------------------------------------------

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
function beam_mass(stations::Vector{Station}, props::Vector{SectionProperties},
                   mat_box::MaterialProperties, mat_reinforcer::MaterialProperties)
    n   = length(stations)
    n == 0 && return 0.0, Float64[]

    # Per-station cross-sectional areas
    width  = [p.width       for p in props]
    height = [p.height      for p in props]
    t_box  = [s.thickness   for s in stations]
    t_ref  = [s.reinforcer  for s in stations]
    spans  = [s.span        for s in stations]

    # Element areas: two flanges + two webs (box) + two reinforcer booms
    web_length = [p.web_length for p in props]
    A_cf  = t_box .* width          # compression flange
    A_tf  = t_box .* width          # tension flange
    A_fw  = t_box .* web_length     # front web
    A_bw  = t_box .* web_length     # back web
    A_ref = 2 .* t_ref.^2           # two reinforcer squares

    A_box = A_cf .+ A_tf .+ A_fw .+ A_bw
    A_tot = A_box .+ A_ref

    # Volume by trapezoidal integration over span
    ds  = diff(spans)
    V_box = 0.5 .* ds .* (A_box[1:end-1] .+ A_box[2:end])
    V_ref = 0.5 .* ds .* (A_ref[1:end-1] .+ A_ref[2:end])

    dm   = V_box .* mat_box.dens .+ V_ref .* mat_reinforcer.dens
    mass = sum(dm)
    return mass, dm
end

# ---------------------------------------------------------------------------
# Report 45 – wing flexibility factor
# ---------------------------------------------------------------------------

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
function report_45(stations::Vector{Station}, props::Vector{SectionProperties},
                   wing::WingGeometry, mat::MaterialProperties)
    M = 1.0          # unit torque (N·m)
    M_imp = M * 0.73756  # N·m → ft·lb

    hasfield(typeof(wing), :V_d) || return 0.0   # graceful degradation
    Vd = getfield(wing, :V_d) * 2.23694          # m/s → mph

    # Unit-torque torsion angle
    Mt_unit = fill(M, length(stations))
    phi = wing_torsion(stations, props, Mt_unit, mat)

    spans = [s.span for s in stations]
    chords = [s.chord for s in stations]

    # Aileron region: start at the span transition (second break)
    s_ail = wing.span_transitions[2]
    i_ail = findfirst(sp -> sp >= s_ail, spans)
    isnothing(i_ail) && (i_ail = 1)

    s_seg  = spans[i_ail:end]
    C_ft   = chords[i_ail:end] .* 3.28084
    theta  = phi[i_ail:end] ./ M_imp

    n = length(s_seg)
    Δs = diff(vcat(s_seg, [wing.span_transitions[end]])) .* 3.28084

    F = sum(theta[1:n] .* C_ft.^2 .* Δs)
    criterion = 200.0 / Vd^2
    return F / criterion
end
