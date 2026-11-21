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
        throw(DimensionMismatch("aerodynamic load arrays must have matching lengths"))
    !isempty(section_props) || throw(ArgumentError("section_props must not be empty"))
    wing.full_span > 0 && wing.constraint_distance > 0 ||
        throw(ArgumentError("span and support distance must be positive"))
    positive = findall(>=(0), load.x)
    if length(positive) != length(section_props)
        positive = findall(>(0), load.x)
    end
    length(positive) == length(section_props) || throw(
        DimensionMismatch(
            "one section property is required for every right-wing load panel",
        ),
    )
    indices = sort(positive; by = i -> load.x[i])
    panel_x = Float64.(load.x[indices])
    halfspan = wing.full_span / 2
    all(z -> 0 <= z <= halfspan, panel_x) ||
        throw(ArgumentError("load positions must lie within the semi-span"))
    issorted(panel_x) ||
        throw(ArgumentError("load positions must be sortable along the semi-span"))
    support = wing.constraint_distance / 2
    support <= halfspan ||
        throw(ArgumentError("strut support must lie within the semi-span"))
    # Interpolate center-of-panel loads and integrate over the half-wing.
    q = 1.5 * 1.15 .* Float64.(load.dLdb[indices])
    m = 1.5 * 1.15 .* Float64.(load.dMadb[indices])
    x = sort(unique(vcat(0.0, panel_x, halfspan)))
    function interp_panel(values, z)
        length(panel_x) == 1 && return values[1]
        z <= panel_x[1] && return values[1]
        z >= panel_x[end] && return values[end]
        j = searchsortedlast(panel_x, z)
        return values[j] +
               (z - panel_x[j]) * (values[j+1] - values[j]) / (panel_x[j+1] - panel_x[j])
    end
    qx = [interp_panel(q, z) for z in x]
    mx = [interp_panel(m, z) for z in x]
    total_lift = trapz(x, qx)
    total_moment = trapz(x, x .* qx)
    reaction_support = total_moment / support
    reaction_root = total_lift - reaction_support
    load_at(z, values) = begin
        z >= halfspan && return 0.0
        zz = sort(unique(vcat(z, x[x .> z], halfspan)))
        trapz(zz, [interp_panel(values, u) for u in zz])
    end
    load_moment(z, values) = begin
        z >= halfspan && return 0.0
        zz = sort(unique(vcat(z, x[x .> z], halfspan)))
        trapz(zz, [interp_panel(values, u) * (u - z) for u in zz])
    end
    xs = panel_x
    v_right = [load_at(z, q) - (z <= support ? reaction_support : 0.0) for z in xs]
    m_right = [
        load_moment(z, q) - (z <= support ? reaction_support * (support - z) : 0.0) for
        z in xs
    ]
    centroid_x = [section_props[findfirst(==(z), panel_x)].centroid[1] for z in xs]
    torsion_dist = m .+ centroid_x .* q
    mt_right = [load_at(z, torsion_dist) for z in xs]
    axial = [
        z <= support ? reaction_support * tan(deg2rad(wing.constraint_angle)) : 0.0 for
        z in xs
    ]
    right = HalfLoads(xs, v_right, m_right, mt_right, axial, argmin(abs.(xs .- support)))
    left = HalfLoads(
        -reverse(xs),
        -reverse(v_right),
        reverse(m_right),
        -reverse(mt_right),
        reverse(axial),
        length(xs) - right.i + 1,
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
    constraint_type == "ss" && return internal_loads_SS(load, section_props, wing)
    constraint_type == "fx" && return internal_loads_FX(load, section_props, wing)
    throw(ArgumentError("unsupported constraint type \"$constraint_type\"; choose \"ss\""))
end
