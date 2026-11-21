"""
Internal structural loads.
Translates `forces/internal_loads_SS.m`.

`internal_loads_FX.m` (biplane/fixed support model) is stubbed; it contained
syntax errors and was work-in-progress in the original Matlab codebase.
"""

# ---------------------------------------------------------------------------
# Simply-supported wing with optional strut support  (internal_loads_SS)
# ---------------------------------------------------------------------------

"""
    internal_loads_SS(load, section_props, wing) -> InternalLoads

Compute internal loads (shear V, bending moment M, torsion Mt, axial F)
for a symmetric wing with two strut supports.

Arguments:
- `load`          – AeroLoad (full-span symmetric distribution)
- `section_props` – Vector{SectionProperties} for the RIGHT semi-span panels
- `wing`          – WingGeometry

Returns InternalLoads with `.L` (left) and `.R` (right) HalfLoads.
"""
function internal_loads_SS(
    load::AeroLoad,
    section_props::Vector{SectionProperties},
    wing::WingGeometry,
)
    n = length(load.x)
    n >= 2 && iseven(n) ||
        throw(ArgumentError("aerodynamic loads must contain an even number of panels"))
    length(load.dLdb) == n && length(load.dMadb) == n ||
        throw(DimensionMismatch("aerodynamic load vectors must match load.x"))
    length(section_props) == n ÷ 2 ||
        throw(DimensionMismatch("one section property is required per right-half panel"))
    all(isfinite, load.x) && all(isfinite, load.dLdb) && all(isfinite, load.dMadb) ||
        throw(ArgumentError("aerodynamic load data must be finite"))
    all(diff(load.x) .> 0) ||
        throw(ArgumentError("aerodynamic panel positions must be strictly increasing"))
    all(isapprox(load.x[i], -load.x[end+1-i]; atol = 1e-8, rtol = 1e-8) for i = 1:(n÷2)) ||
        throw(ArgumentError("aerodynamic panel positions must be symmetric about the root"))
    all(prop -> !isempty(prop.centroid) && isfinite(prop.centroid[1]), section_props) ||
        throw(ArgumentError("section centroids must contain finite x coordinates"))

    b = wing.full_span
    d = wing.constraint_distance
    isfinite(b) && b > 0 || throw(ArgumentError("full span must be positive and finite"))
    0 < d < b || throw(ArgumentError("support spacing must be between zero and full span"))
    isfinite(wing.constraint_angle) || throw(ArgumentError("support angle must be finite"))
    load.x[1] > -b / 2 && load.x[end] < b / 2 ||
        throw(ArgumentError("panel centers must lie strictly inside the wing tips"))

    # MATLAB's load arrays are panel-centred; the zero-valued tips complete
    # the trapezoidal integration intervals.
    s = vcat(-b / 2, load.x, b / 2)
    s0 = s .+ b / 2
    factor = 1.5 * 1.15
    w = vcat(0.0, load.dLdb .* factor, 0.0)
    wt = vcat(0.0, load.dMadb .* factor, 0.0)

    # Section centroids are measured from the leading edge. Positive
    # eccentricity is aft of the quarter-chord, as in the MATLAB model.
    chord =
        x -> begin
            y = abs(x)
            breaks = wing.span_transitions
            chords = wing.chord_transitions
            length(breaks) == length(chords) && length(breaks) >= 2 ||
                throw(ArgumentError("wing chord and span transitions must match"))
            all(isfinite, breaks) && all(isfinite, chords) && all(>(0), chords) ||
                throw(
                    ArgumentError(
                        "wing transitions and chords must be finite and positive",
                    ),
                )
            all(diff(breaks) .> 0) || throw(
                ArgumentError("wing span transitions must be strictly increasing"),
            )
            breaks[1] <= y <= breaks[end] ||
                throw(ArgumentError("panel station lies outside wing span transitions"))
            i = min(searchsortedlast(breaks, y), length(breaks) - 1)
            α = (y - breaks[i]) / (breaks[i+1] - breaks[i])
            (1 - α) * chords[i] + α * chords[i+1]
        end
    right_indices = (n÷2+1):n
    e_right = [
        prop.centroid[1] - chord(load.x[j]) / 4 for
        (prop, j) in zip(section_props, right_indices)
    ]
    e = vcat(reverse(e_right), e_right)
    arm_load = vcat(0.0, e .* load.dLdb .* factor .+ load.dMadb .* factor, 0.0)

    # The source model locates reactions at the nearest span stations. The
    # support angle is measured from vertical (0° is a vertical strut).
    i1 = argmin(abs.(s0 .- (b - d) / 2))
    i2 = argmin(abs.(s0 .- (b + d) / 2))
    1 < i1 < i2 < length(s) ||
        throw(ArgumentError("support stations must be distinct and inside the wing"))
    θ = deg2rad(wing.constraint_angle)
    cosθ = cos(θ)
    abs(cosθ) > eps(Float64) || throw(
        ArgumentError("horizontal supports cannot balance vertical aerodynamic loads"),
    )

    int_w = cumtrapz(s0, w)
    int_ws = cumtrapz(s0, w .* s0)
    reactions = [
        cosθ cosθ
        cosθ * s0[i1] cosθ * s0[i2]
    ] \ [int_w[end], int_ws[end]]
    Ra, Rb = reactions
    Ra_v, Rb_v = Ra * cosθ, Rb * cosθ

    V1 = int_w
    M1 = s0 .* V1 .- int_ws
    V2 = int_w .- Ra_v
    M2 = s0 .* (V2 .+ Ra_v) .- int_ws .- Ra_v .* (s0 .- s0[i1])
    V3 = int_w .- Ra_v .- Rb_v
    M3 =
        s0 .* (V3 .+ Ra_v .+ Rb_v) .- int_ws .- Ra_v .* (s0 .- s0[i1]) .-
        Rb_v .* (s0 .- s0[i2])

    arm1 = cumtrapz(s0[1:i1], arm_load[1:i1])
    arm2 = cumtrapz(s0[i1:i2], arm_load[i1:i2])
    arm3 = cumtrapz(s0[i2:end], arm_load[i2:end])
    torque2 = trapz(s0[i1:i2], arm_load[i1:i2])
    x_torque =
        isapprox(torque2, 0.0; atol = eps(Float64)) ? (s0[i1] + s0[i2]) / 2 :
        trapz(s0[i1:i2], arm_load[i1:i2] .* s0[i1:i2]) / torque2
    Mt1 = arm1
    Mt2 = -arm2 .+ torque2 .* (s0[i2] - x_torque) / d
    Mt3 = arm3 .- arm3[end]

    axial = range(Ra * sin(θ), Rb * sin(θ); length = i2 - i1 + 1)
    join_torque(a, b) = begin
        magnitude = max(abs(a), abs(b))
        direction = iszero(a) ? (iszero(b) ? 1.0 : sign(b)) : sign(a)
        direction * magnitude
    end
    V = vcat(
        V1[1:(i1-1)],
        abs(V1[i1] - V2[i1]),
        V2[(i1+1):(i2-1)],
        -abs(V2[i2] - V3[i2]),
        V3[(i2+1):end],
    )
    M = vcat(
        M1[1:(i1-1)],
        max(M1[i1], M2[i1]),
        M2[(i1+1):(i2-1)],
        max(M2[i2], M3[i2]),
        M3[(i2+1):end],
    )
    Mt = vcat(
        Mt1[1:(end-1)],
        join_torque(Mt1[end], Mt2[1]),
        Mt2[2:(end-1)],
        join_torque(Mt2[end], Mt3[1]),
        Mt3[2:end],
    )
    F = vcat(zeros(i1 - 1), collect(axial), zeros(length(s) - i2))

    # Drop the artificial zero-load tip stations, then orient each half from
    # root to tip. Left shear and torsion use the opposite local sign.
    s = s[2:(end-1)]
    V, M, Mt, F = V[2:(end-1)], M[2:(end-1)], Mt[2:(end-1)], F[2:(end-1)]
    mid = length(s) ÷ 2
    left = HalfLoads(
        reverse(-s[1:mid]),
        reverse(-V[1:mid]),
        reverse(M[1:mid]),
        reverse(-Mt[1:mid]),
        reverse(F[1:mid]),
        mid - i1 + 1,
    )
    right = HalfLoads(
        s[(mid+1):end],
        V[(mid+1):end],
        M[(mid+1):end],
        Mt[(mid+1):end],
        F[(mid+1):end],
        i2 - (length(s0) ÷ 2 + 1),
    )
    return InternalLoads(left, right)
end

# ---------------------------------------------------------------------------
# Fixed-support (FX) stub
# ---------------------------------------------------------------------------

"""
    internal_loads_FX(load, section_props, wing) -> InternalLoads

**Not yet implemented.** The original `internal_loads_FX.m` contains syntax
errors and was a work-in-progress biplane model.
"""
function internal_loads_FX(::AeroLoad, ::Vector{SectionProperties}, ::WingGeometry)
    error(
        "internal_loads_FX is not yet implemented. " *
        "Use constraint_type = \"ss\" for simply-supported analysis.",
    )
end

# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

"""
    compute_internal_loads(load, section_props, wing, constraint_type) -> InternalLoads
"""
function compute_internal_loads(
    load::AeroLoad,
    section_props::Vector{SectionProperties},
    wing::WingGeometry,
    constraint_type::String,
)
    constraint = lowercase(constraint_type)
    constraint == "ss" && return internal_loads_SS(load, section_props, wing)
    constraint == "fx" && return internal_loads_FX(load, section_props, wing)
    throw(ArgumentError("unknown constraint type: $constraint_type"))
end
