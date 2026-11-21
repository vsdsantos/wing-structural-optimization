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
    geometry = foil_properties(station.chord, 0.0, station.limits, [0.0, 0.0], foil)
    upper_front, upper_back = geometry.upper[1, :], geometry.upper[end, :]
    lower_front, lower_back = geometry.lower[1, :], geometry.lower[end, :]

    A = upper_front[2] < upper_back[2] ?
        upper_front : [upper_front[1], upper_back[2]]
    B = upper_front[2] < upper_back[2] ?
        [upper_back[1], upper_front[2]] : upper_back
    C = lower_front[2] < lower_back[2] ?
        [lower_front[1], lower_back[2]] : lower_front
    D = lower_front[2] < lower_back[2] ?
        lower_back : [lower_back[1], lower_front[2]]

    width = D[1] - C[1]
    height = A[2] - C[2]
    thickness = station.thickness
    web_length = height - 2thickness
    width > thickness > 0 ||
        throw(
            DomainError(
                (width, thickness),
                "box width must exceed positive wall thickness",
            ),
        )
    web_length > 0 || throw(DomainError(web_length, "wall thickness leaves no clear web"))

    flange_area = thickness * width
    web_area = thickness * web_length
    compression_centroid = A .+ [width, -thickness] ./ 2
    tension_centroid = C .+ [width, thickness] ./ 2
    front_web_centroid = C .+ [thickness / 2, thickness + web_length / 2]
    back_web_centroid = D .+ [-thickness / 2, thickness + web_length / 2]
    element_centroids = [
        compression_centroid'
        tension_centroid'
        front_web_centroid'
        back_web_centroid'
    ]
    element_areas = [flange_area, flange_area, web_area, web_area]
    centroid = multielement_centroid(element_centroids, element_areas)

    web_inertia = [thickness * web_length^3 / 12 0.0; 0.0 thickness^3 * web_length / 12]
    flange_inertia = [width * thickness^3 / 12 0.0; 0.0 width^3 * thickness / 12]
    inertia =
        translate_inertia_tensor(
            flange_inertia, compression_centroid - centroid, flange_area
        ) +
        translate_inertia_tensor(
            flange_inertia, tension_centroid - centroid, flange_area
        ) +
        translate_inertia_tensor(web_inertia, front_web_centroid - centroid, web_area) +
        translate_inertia_tensor(web_inertia, back_web_centroid - centroid, web_area)
    inertia[1, 2] = inertia[2, 1] = 0.0
    principal_inertia, theta = principal_inertia_tensor(inertia)
    vertices = [A'; B'; C'; D']
    vertices_mid = vertices .+ thickness .* [1 -1; -1 -1; 1 1; -1 1] ./ 2
    area_mid = (width - thickness) * (height - thickness)
    area = sum(element_areas)

    return SectionProperties(
        inertia, principal_inertia, principal_inertia[1, 1], principal_inertia[2, 2],
        theta, centroid, area_mid, area, width, height, web_length, vertices,
        vertices_mid, compression_centroid, tension_centroid,
    )
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
    throw(
        ArgumentError(
            "D-section properties are not supported by the current SectionProperties model",
        ),
    )
end

# ---------------------------------------------------------------------------
# 4. O-section (circular)
# ---------------------------------------------------------------------------

"""
    O_section(station, foil; offset=0.0) -> SectionProperties

Circular tube inscribed in the airfoil section at the spar limits.
"""
function O_section(station::Station, foil::FoilData; offset::Float64 = 0.0)
    throw(
        ArgumentError(
            "O-section properties are not supported by the current SectionProperties model",
        ),
    )
end

# ---------------------------------------------------------------------------
# 5. Foil-shaped section
# ---------------------------------------------------------------------------

"""
    foil_shaped_section(station, foil) -> SectionProperties

Full foil-profile section: curved upper/lower flanges + front and back webs.
"""
function foil_shaped_section(station::Station, foil::FoilData)
    throw(
        ArgumentError(
            "foil-shaped section properties are not supported",
        ),
    )
end

# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

"""
    compute_section(station, foil, section_type) -> SectionProperties

Dispatch to the correct section-computation function based on `section_type`.
"""
function compute_section(station::Station, foil::FoilData, section_type::String)
    section_type in ("box", "reinforced_box") ||
        throw(ArgumentError("unsupported section type: $section_type"))
    return reinforced_box_section(station, foil)
end
