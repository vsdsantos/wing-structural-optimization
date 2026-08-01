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
    error("Not implemented")
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
    error("Not implemented")
end

# ---------------------------------------------------------------------------
# 4. O-section (circular)
# ---------------------------------------------------------------------------

"""
    O_section(station, foil; offset=0.0) -> SectionProperties

Circular tube inscribed in the airfoil section at the spar limits.
"""
function O_section(station::Station, foil::FoilData; offset::Float64 = 0.0)
    error("Not implemented")
end

# ---------------------------------------------------------------------------
# 5. Foil-shaped section
# ---------------------------------------------------------------------------

"""
    foil_shaped_section(station, foil) -> SectionProperties

Full foil-profile section: curved upper/lower flanges + front and back webs.
"""
function foil_shaped_section(station::Station, foil::FoilData)
    error("Not implemented")
end

# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

"""
    compute_section(station, foil, section_type) -> SectionProperties

Dispatch to the correct section-computation function based on `section_type`.
"""
function compute_section(station::Station, foil::FoilData, section_type::String)
    error("Not implemented")
end
