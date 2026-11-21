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
for a symmetric wing with a strut (simply-supported) constraint.
Support reactions use the nearest sampled support stations and root-referenced
equilibrium; this removes the mixed-origin reaction balance in the MATLAB
script and can change results slightly when a support falls between stations.

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
    length(load.x) == length(load.dLdb) == length(load.dMadb) ||
        throw(
            DimensionMismatch(
                "aerodynamic positions, lift, and pitching moment must have equal lengths",
            ),
        )
    !isempty(load.x) ||
        throw(ArgumentError("at least one aerodynamic load station is required"))
    issorted(load.x) ||
        throw(ArgumentError("aerodynamic stations must be ordered from left to right"))
    b = wing.full_span
    all(abs.(load.x) .<= b / 2) ||
        throw(DomainError(load.x, "aerodynamic stations must lie within the full span"))
    0 < wing.constraint_distance < b ||
        throw(
            DomainError(
                wing.constraint_distance,
                "support spacing must be between zero and the full span",
            ),
        )
    length(section_props) * 2 == length(load.x) ||
        throw(
            DimensionMismatch(
                "section_props must describe the right half of the aerodynamic stations",
            ),
        )

    s = [-b / 2; load.x; b / 2]
    s0 = s .+ b / 2
    n = length(s)
    w = [0.0; 1.5 * 1.15 .* load.dLdb; 0.0]
    wt = [0.0; 1.5 * 1.15 .* load.dMadb; 0.0]
    centroids = [prop.centroid[1] for prop in section_props]
    xl = [centroids[end]; reverse(centroids); centroids; centroids[end]]

    nominal_x1 = (b - wing.constraint_distance) / 2
    nominal_x2 = (b + wing.constraint_distance) / 2
    i1 = argmin(abs.(s0 .- nominal_x1))
    i2 = argmin(abs.(s0 .- nominal_x2))
    1 < i1 < i2 < n ||
        throw(ArgumentError("support stations must lie inside the aerodynamic span"))
    x1, x2 = s0[i1], s0[i2]
    support_span = x2 - x1
    support_span > 0 || throw(ArgumentError("support stations must be distinct"))

    int_w = cumtrapz(s0, w)
    int_ws = cumtrapz(s0, w .* s0)
    int_wt1 = cumtrapz(s0[1:i1], (wt .+ xl .* w)[1:i1])
    int_wt2 = cumtrapz(s0[i1:i2], (wt .+ xl .* w)[i1:i2])
    int_wt3 = cumtrapz(s0[i2:end], (wt .+ xl .* w)[i2:end])

    # Resolve the two simple supports from vertical-force and root-moment equilibrium.
    theta = deg2rad(wing.constraint_angle)
    vertical_factor = cos(theta)
    abs(vertical_factor) > eps(Float64) ||
        throw(DomainError(wing.constraint_angle, "support angle cannot be horizontal"))
    reaction_matrix = [1.0 1.0; x1 x2]
    Ra, Rb = reaction_matrix \ [int_w[end], int_ws[end]]
    axial_reaction_a, axial_reaction_b = Ra * tan(theta), Rb * tan(theta)

    V1 = int_w
    M1 = s0 .* V1 .- int_ws
    torque1 = int_wt1

    V2 = int_w .- Ra
    M2 = s0 .* (V2 .+ Ra) .- int_ws .- Ra .* (s0 .- x1)
    wt_section = wt[i1:i2]
    torque_denominator = trapz(s0[i1:i2], wt_section)
    torque_centroid = iszero(torque_denominator) ? (x1 + x2) / 2 :
                      trapz(s0[i1:i2], wt_section .* s0[i1:i2]) / torque_denominator
    torque2 = -int_wt2 .+ int_wt2[end] .* (x2 - torque_centroid) ./ support_span

    V3 = int_w .- Ra .- Rb
    M3 = s0 .* (V3 .+ Ra .+ Rb) .- int_ws .- Ra .* (s0 .- x1) .- Rb .* (s0 .- x2)
    torque3 = int_wt3 .- int_wt3[end]
    N2 = range(axial_reaction_a, axial_reaction_b; length = i2 - i1 + 1)

    loads_V = [
        V1[1:(i1 - 1)];
        abs(V1[i1] - V2[i1]);
        V2[(i1 + 1):(i2 - 1)];
        -abs(V2[i2] - V3[i2]);
        V3[(i2 + 1):end]
    ]
    loads_M = [
        M1[1:(i1 - 1)];
        max(M1[i1], M2[i1]);
        M2[(i1 + 1):(i2 - 1)];
        max(M2[i2], M3[i2]);
        M3[(i2 + 1):end]
    ]
    blend(a, b) = iszero(a) ? sign(b) * abs(b) : sign(a) * max(abs(a), abs(b))
    loads_Mt = [
        torque1[1:(end - 1)];
        blend(torque1[end], torque2[1]);
        torque2[2:(end - 1)];
        blend(torque3[1], torque2[end]);
        torque3[2:end]
    ]
    loads_F = [zeros(i1 - 1); collect(N2); zeros(n - i2)]
    lengths = (length(loads_V), length(loads_M), length(loads_Mt), length(loads_F))
    all(==(n), lengths) ||
        throw(
            DimensionMismatch(
                "piecewise internal-load assembly produced inconsistent station counts",
            ),
        )

    span_stations = s[2:(end - 1)]
    loads_V = loads_V[2:(end - 1)]
    loads_M = loads_M[2:(end - 1)]
    loads_Mt = loads_Mt[2:(end - 1)]
    loads_F = loads_F[2:(end - 1)]
    midpoint = (n - 2) ÷ 2
    iseven(n - 2) ||
        throw(
            ArgumentError(
                "symmetric loading requires an even number of non-tip stations",
            ),
        )
    right_range = (midpoint + 1):(n - 2)
    left_range = midpoint:-1:1
    right_index = clamp(i2 - midpoint - 1, 1, length(right_range))
    left_index = clamp(midpoint - i1 + 1, 1, midpoint)
    return InternalLoads(
        HalfLoads(
            -span_stations[left_range], -loads_V[left_range], loads_M[left_range],
            -loads_Mt[left_range], loads_F[left_range], left_index,
        ),
        HalfLoads(
            span_stations[right_range], loads_V[right_range], loads_M[right_range],
            loads_Mt[right_range], loads_F[right_range], right_index,
        ),
    )
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
    if constraint_type == "ss"
        return internal_loads_SS(load, section_props, wing)
    elseif constraint_type == "fx"
        return internal_loads_FX(load, section_props, wing)
    end
    throw(ArgumentError("constraint_type must be \"ss\" or \"fx\""))
end
