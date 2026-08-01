"""
Visualization functions.
Translates `plots/`: plot_esforcos.m, plot_longarina.m,
plot_prop_geo.m, plot_sec_retangular.m (reinforced box), plot_sec_perfilada.m.

All functions return the Plots.jl plot object so the caller can save/display.
"""

"""
    plot_internal_loads(loads_half; title_prefix="Right wing") -> Plots.Plot

Plot shear force, bending moment, torsional moment, and axial force
diagrams for one wing half.
"""
function plot_internal_loads(hl::HalfLoads; title_prefix::String = "Right wing")
    error("Not implemented")
end

"""
    plot_stiffness(stations, props, mat) -> Plots.Plot

Plot bending stiffness EI and torsional stiffness GJ vs semi-span.
"""
function plot_stiffness(
    stations::Vector{Station},
    props::Vector{SectionProperties},
    mat::MaterialProperties,
)
    error("Not implemented")
end

"""
    plot_box_section(prop; title="Box section") -> Plots.Plot

Draw the cross-section outline of a reinforced box section.
"""
function plot_box_section(prop::SectionProperties; title::String = "Box section")
    error("Not implemented")
end

"""
    plot_spar_topview(stations, wing) -> Plots.Plot

Top-view sketch of the wing planform with the spar overlay.
"""
function plot_spar_topview(stations::Vector{Station}, wing::WingGeometry)
    error("Not implemented")
end
