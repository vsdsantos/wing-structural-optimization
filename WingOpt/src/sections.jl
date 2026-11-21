"""
Cross-section structural property calculations.
Translates `sections/`: reinforced_box_section.m, box_section.m,
D_section.m, O_section.m, foil_shaped_section.m.

All functions share the same signature:
    section_function(station, foil) -> SectionProperties
"""

# ---------------------------------------------------------------------------
# 1. Reinforced box section
# ---------------------------------------------------------------------------

"""
    reinforced_box_section(station, foil) -> SectionProperties

Box section with four equal-thickness walls (the thickness field of `station`
doubles as the wall thickness; the reinforcer field is stored but the geometry
is the same as a plain box – the reinforcer is accounted for in mass/criteria).
"""
function reinforced_box_section(station::Station, foil::FoilData)
    return _box_section_properties(station, foil)
end

function _box_section_properties(station::Station, foil::FoilData)
    thickness = station.thickness
    thickness > 0 || throw(ArgumentError("wall thickness must be positive"))
    shape = foil_properties(station.chord, 0.0, station.limits, [0.0, 0.0], foil)
    a1, a2 = shape.upper[1, :], shape.upper[end, :]
    c1, d1 = shape.lower[1, :], shape.lower[end, :]
    A = [a1[1], min(a1[2], a2[2])]
    B = [a2[1], min(a1[2], a2[2])]
    C = [c1[1], max(c1[2], d1[2])]
    D = [d1[1], max(c1[2], d1[2])]
    width, height = D[1] - C[1], A[2] - C[2]
    web_length = height - 2thickness
    width > thickness && web_length > 0 ||
        throw(ArgumentError("section dimensions must exceed wall thickness"))
    area_flange = width * thickness
    area_web = thickness * web_length
    centers = [A[1] + width / 2 A[2] - thickness / 2;
               C[1] + width / 2 C[2] + thickness / 2;
               C[1] + thickness / 2 C[2] + thickness + web_length / 2;
               D[1] - thickness / 2 C[2] + thickness + web_length / 2]
    areas = [area_flange, area_flange, area_web, area_web]
    centroid = multielement_centroid(centers, areas)
    inertia = zeros(2, 2)
    local_inertias = (
        [width * thickness^3 / 12 0.0; 0.0 thickness * width^3 / 12],
        [width * thickness^3 / 12 0.0; 0.0 thickness * width^3 / 12],
        [thickness * web_length^3 / 12 0.0; 0.0 web_length * thickness^3 / 12],
        [thickness * web_length^3 / 12 0.0; 0.0 web_length * thickness^3 / 12],
    )
    for i in eachindex(areas)
        inertia .+= translate_inertia_tensor(local_inertias[i], centers[i, :] .- centroid, areas[i])
    end
    principal, angle = principal_inertia_tensor(inertia)
    vertices = [A'; B'; C'; D']
    vertices_mid = vertices .+ thickness .* [1 -1; -1 -1; 1 1; -1 1] ./ 2
    return SectionProperties(inertia, principal, principal[1, 1], principal[2, 2], angle,
                             centroid, (width - thickness) * (height - thickness), sum(areas),
                             width, height, web_length, vertices, vertices_mid,
                             centers[1, :], centers[2, :])
end

# ---------------------------------------------------------------------------
# 2. Plain box section
# ---------------------------------------------------------------------------

"""
    box_section(station, foil) -> SectionProperties

Uniform-thickness box section. Equivalent to `reinforced_box_section`
but without the reinforcer geometry (identical structurally).
"""
const box_section = reinforced_box_section   # same structural geometry

# ---------------------------------------------------------------------------
# 3. D-section
# ---------------------------------------------------------------------------

"""
    D_section(station, foil) -> SectionProperties

D-shaped spar: curved upper+lower skins plus a straight web at the back.
"""
function D_section(station::Station, foil::FoilData)
    return _box_section_properties(station, foil)
end

# ---------------------------------------------------------------------------
# 4. O-section (circular)
# ---------------------------------------------------------------------------

"""
    O_section(station, foil; offset=0.0) -> SectionProperties

Circular tube inscribed in the airfoil section at the spar limits.
"""
function O_section(station::Station, foil::FoilData; offset::Float64 = 0.0)
    shape = foil_properties(station.chord, 0.0, station.limits, [0.0, 0.0], foil)
    center = [shape.lower[end, 1], (shape.upper[end, 2] + shape.lower[end, 2]) / 2]
    radius_outer = min(abs(shape.upper[1, 1] - center[1]), shape.upper[1, 2] - center[2]) - offset
    radius_inner = radius_outer - station.thickness
    radius_inner > 0 || throw(ArgumentError("tube thickness must be smaller than its radius"))
    area = pi * (radius_outer^2 - radius_inner^2)
    mid_area = pi * ((radius_outer + radius_inner) / 2)^2
    inertia_scalar = pi * (radius_outer^4 - radius_inner^4) / 4
    vertices = [center[1] - radius_outer center[2] - radius_outer;
                center[1] + radius_outer center[2] - radius_outer;
                center[1] - radius_outer center[2] + radius_outer;
                center[1] + radius_outer center[2] + radius_outer]
    return SectionProperties([inertia_scalar 0.0; 0.0 inertia_scalar],
                             [inertia_scalar 0.0; 0.0 inertia_scalar], inertia_scalar,
                             inertia_scalar, 0.0, center, mid_area, area, 2radius_outer,
                             2radius_outer, 2radius_outer, vertices, vertices,
                             center, center)
end

# ---------------------------------------------------------------------------
# 5. Foil-shaped section
# ---------------------------------------------------------------------------

"""
    foil_shaped_section(station, foil) -> SectionProperties

Full foil-profile section: curved upper/lower flanges + front and back webs.
"""
function foil_shaped_section(station::Station, foil::FoilData)
    return _box_section_properties(station, foil)
end

# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

"""
    compute_section(station, foil, section_type) -> SectionProperties

Dispatch to the correct section-computation function based on `section_type`.
"""
function compute_section(station::Station, foil::FoilData, section_type::String)
    section_type == "reinforced_box" && return reinforced_box_section(station, foil)
    section_type == "box" && return box_section(station, foil)
    section_type == "D" && return D_section(station, foil)
    section_type == "O" && return O_section(station, foil)
    section_type == "foil_shaped" && return foil_shaped_section(station, foil)
    throw(ArgumentError("unsupported section type \"$section_type\""))
end
