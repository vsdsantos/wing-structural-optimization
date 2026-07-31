"""
Cross-section structural property calculations.
Translates `sections/`: reinforced_box_section.m, box_section.m,
D_section.m, O_section.m, foil_shaped_section.m.

All functions share the same signature:
    section_function(station, foil) -> SectionProperties

Notes on known issues fixed during migration:
- box_section.m used old Portuguese helper names (centroide_multielemento,
  translada_eixo, inercia_principal) – unified to the English API here.
- foil_shaped_section (typo "foild_shaped_section" in fit_function.m) corrected.
"""

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""Find the axis-aligned bounding box of the foil cross-section."""
function _inscribed_rect(fg::FoilGeometry)
    a1_c = fg.upper[1, :];   a2_c = fg.upper[end, :]
    a1_t = fg.lower[1, :];   a2_t = fg.lower[end, :]

    # Top-left (A) and top-right (B) corners of the compression flange
    if a1_c[2] < a2_c[2]
        A = copy(a1_c);              B = [a2_c[1], a1_c[2]]
    else
        A = [a1_c[1], a2_c[2]];     B = copy(a2_c)
    end

    # Bottom-left (C) and bottom-right (D) corners of the tension flange
    if a1_t[2] < a2_t[2]
        C = [a1_t[1], a2_t[2]];     D = copy(a2_t)
    else
        C = copy(a1_t);              D = [a2_t[1], a1_t[2]]
    end

    return A, B, C, D
end

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
    chord = station.chord;  t = station.thickness
    fg = foil_properties(chord, 0.0, station.limits, [0.0, 0.0], foil)
    A, B, C, D = _inscribed_rect(fg)

    width      = D[1] - C[1]
    height     = A[2] - C[2]
    web_length = height - 2*t

    # Element areas
    a_cf = t * width          # compression flange
    a_tf = t * width          # tension flange
    a_w  = t * web_length     # one web (×2)

    # Element centroids
    c_cf = A .+ [width, -t] ./ 2
    c_tf = C .+ [width,  t] ./ 2
    c_fw = C .+ [t/2,   t + web_length/2]
    c_bw = D .+ [-t/2,  t + web_length/2]

    areas = [a_cf, a_tf, a_w, a_w]
    cents = vcat(c_cf', c_tf', c_fw', c_bw')
    centroid = multielement_centroid(cents, areas)

    # Second moments about own centroidal axes (rectangles)
    I_web = [t*web_length^3/12  0.0; 0.0  t^3*web_length/12]
    I_cf  = [width*t^3/12       0.0; 0.0  width^3*t/12      ]
    I_tf  = [width*t^3/12       0.0; 0.0  width^3*t/12      ]

    # Parallel-axis theorem → section centroid
    I_cf = translate_inertia_tensor(I_cf, c_cf .- centroid, a_cf)
    I_tf = translate_inertia_tensor(I_tf, c_tf .- centroid, a_tf)
    I_fw = translate_inertia_tensor(I_web, c_fw .- centroid, a_w)
    I_bw = translate_inertia_tensor(I_web, c_bw .- centroid, a_w)

    I = I_cf .+ I_tf .+ I_fw .+ I_bw
    I[1, 2] = 0.0;  I[2, 1] = 0.0    # symmetric section → Ixy = 0

    I_max, theta = principal_inertia_tensor(I)

    area_mid = (width - t) * (height - t)
    vert     = vcat(A', B', C', D')
    vert_mid = vert .+ t .* [1 -1; -1 -1; 1 1; -1 1] ./ 2

    return SectionProperties(I, I_max, I_max[1,1], I_max[2,2], theta,
                              centroid, area_mid, sum(areas),
                              width, height, web_length,
                              vert, vert_mid, c_cf, c_tf)
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
    chord = station.chord;  t = station.thickness
    fg = foil_properties(chord, 0.0, station.limits, [0.0, 0.0], foil)

    sup_area = fg.extra_length * t
    sup_cent = section_centroid(fg.upper)

    inf_area = fg.intra_length * t
    inf_cent = section_centroid(fg.lower)

    # Web: use back-web height as a rectangle
    alma_len  = fg.back_length
    alma_area = alma_len * t
    alma_cent = [fg.lower[end, 1], fg.lower[end, 2] + alma_len/2]

    areas    = [sup_area, inf_area, alma_area]
    cents    = vcat(sup_cent', inf_cent', alma_cent')
    centroid = multielement_centroid(cents, areas)

    sup_I = area_inertia_tensor(t, fg.upper .- sup_cent')
    inf_I = area_inertia_tensor(t, fg.lower .- inf_cent')
    alm_I = [t*alma_len^3/12  0.0; 0.0  t^3*alma_len/12]

    sup_I = translate_inertia_tensor(sup_I, sup_cent .- centroid, sup_area)
    inf_I = translate_inertia_tensor(inf_I, inf_cent .- centroid, inf_area)
    alm_I = translate_inertia_tensor(alm_I, alma_cent .- centroid, alma_area)

    I = sup_I .+ inf_I .+ alm_I
    I_max, theta = principal_inertia_tensor(I)

    vert     = zeros(4, 2);  vert_mid = zeros(4, 2)
    return SectionProperties(I, I_max, I_max[1,1], I_max[2,2], theta,
                              centroid, fg.area, sum(areas),
                              fg.back_length, fg.front_length, 0.0,
                              vert, vert_mid, sup_cent, inf_cent)
end

# ---------------------------------------------------------------------------
# 4. O-section (circular)
# ---------------------------------------------------------------------------

"""
    O_section(station, foil; offset=0.0) -> SectionProperties

Circular tube inscribed in the airfoil section at the spar limits.
"""
function O_section(station::Station, foil::FoilData; offset::Float64=0.0)
    chord = station.chord;  t = station.thickness
    fg = foil_properties(chord, 0.0, station.limits, [0.0, 0.0], foil)

    # Centre = mid-point of the back web; radius = minimum distance to foil surface
    t_back   = fg.upper[end, 2] - fg.lower[end, 2]
    C        = [fg.lower[end, 1], fg.lower[end, 2] + t_back/2]
    all_pts  = vcat(fg.upper, fg.lower)
    r_all    = sqrt.(sum((all_pts .- C').^2, dims=2)[:])
    r_max    = minimum(r_all)

    r_ex  = r_max - offset
    r_in  = r_ex  - t
    area  = π * (r_ex^2 - r_in^2) / 2
    Am    = π * ((r_ex - r_in)^2) / 2

    I     = π * (r_ex^4 - r_in^4) .* [1/4 1/2; 1/2 1/4]
    I_max, theta = principal_inertia_tensor(I)

    vert = zeros(4, 2)
    return SectionProperties(I, I_max, I_max[1,1], I_max[2,2], theta,
                              C, Am, area,
                              2*r_ex, 2*r_ex, 0.0,
                              vert, vert, C, C)
end

# ---------------------------------------------------------------------------
# 5. Foil-shaped section
# ---------------------------------------------------------------------------

"""
    foil_shaped_section(station, foil) -> SectionProperties

Full foil-profile section: curved upper/lower flanges + front and back webs.
"""
function foil_shaped_section(station::Station, foil::FoilData)
    chord = station.chord
    t_mc  = station.thickness        # compression-flange thickness
    t_mt  = station.reinforcer       # tension-flange thickness (reuse reinforcer field)
    t_alma = station.thickness       # web thickness (same as box thickness)

    # Flanges with half-thickness offset
    m_f  = foil_properties(chord, 0.0, station.limits, [t_mc/2, t_mt/2], foil)
    mc_area = m_f.extra_length * t_mc
    mc_cent = section_centroid(m_f.upper)
    mt_area = m_f.intra_length * t_mt
    mt_cent = section_centroid(m_f.lower)

    # Webs with full offset
    a_f   = foil_properties(chord, 0.0, station.limits, [t_mc, t_mt], foil)
    a1_area = a_f.front_length * t_alma
    a2_area = a_f.back_length  * t_alma
    a1_cent = [a_f.lower[1,   1], a_f.lower[1,   2] + a_f.front_length/2]
    a2_cent = [a_f.lower[end, 1], a_f.lower[end, 2] + a_f.back_length /2]

    areas    = [mc_area, mt_area, a1_area, a2_area]
    cents    = vcat(mc_cent', mt_cent', a1_cent', a2_cent')
    centroid = multielement_centroid(cents, areas)

    mc_I = area_inertia_tensor(t_mc,  m_f.upper .- mc_cent')
    mt_I = area_inertia_tensor(t_mt,  m_f.lower .- mt_cent')
    a1_I = [t_alma*a_f.front_length^3/12  0.0; 0.0  t_alma^3*a_f.front_length/12]
    a2_I = [t_alma*a_f.back_length^3/12   0.0; 0.0  t_alma^3*a_f.back_length /12]

    mc_I = translate_inertia_tensor(mc_I, mc_cent .- centroid, mc_area)
    mt_I = translate_inertia_tensor(mt_I, mt_cent .- centroid, mt_area)
    a1_I = translate_inertia_tensor(a1_I, a1_cent .- centroid, a1_area)
    a2_I = translate_inertia_tensor(a2_I, a2_cent .- centroid, a2_area)

    I = mc_I .+ mt_I .+ a1_I .+ a2_I
    I_max, theta = principal_inertia_tensor(I)

    foil_lm = foil_properties(chord, 0.0, station.limits, [t_mc/2, t_mt/2], foil)
    vert = zeros(4, 2)
    return SectionProperties(I, I_max, I_max[1,1], I_max[2,2], theta,
                              centroid, foil_lm.area, sum(areas),
                              m_f.back_length, m_f.front_length, 0.0,
                              vert, vert, mc_cent, mt_cent)
end

# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

"""
    compute_section(station, foil, section_type) -> SectionProperties

Dispatch to the correct section-computation function based on `section_type`.
"""
function compute_section(station::Station, foil::FoilData, section_type::String)
    if section_type == "reinforced_box"
        return reinforced_box_section(station, foil)
    elseif section_type == "box"
        return box_section(station, foil)
    elseif section_type == "D"
        return D_section(station, foil)
    elseif section_type == "O"
        return O_section(station, foil)
    elseif section_type == "foil_shaped"
        return foil_shaped_section(station, foil)
    else
        error("Unknown section_type: \"$section_type\"")
    end
end
