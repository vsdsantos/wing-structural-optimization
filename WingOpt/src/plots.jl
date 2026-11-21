"""
Visualization functions.
Translates `plots/`: plot_esforcos.m, plot_longarina.m,
plot_prop_geo.m, plot_sec_retangular.m (reinforced box), plot_sec_perfilada.m.

All functions return the Plots.jl plot object so the caller can save/display.
"""

import Plots

"""
    plot_internal_loads(loads_half; title_prefix="Right wing") -> Plots.Plot

Plot shear force, bending moment, torsional moment, and axial force
diagrams for one wing half.
"""
function plot_internal_loads(hl::HalfLoads; title_prefix::String = "Right wing")
    p = Plots.plot(layout = (4, 1), size = (800, 900), legend = false)
    for (i, (values, label)) in enumerate(((hl.V, "Shear (N)"), (hl.M, "Bending (N·m)"),
                                          (hl.Mt, "Torsion (N·m)"), (hl.F, "Axial (N)")))
        Plots.plot!(p[i], hl.s, values; ylabel = label, title = i == 1 ? title_prefix : "")
    end
    Plots.plot!(p[4]; xlabel = "Span (m)")
    return p
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
    length(stations) == length(props) ||
        throw(DimensionMismatch("stations and properties must match"))
    span = getfield.(stations, :span)
    bending = mat.axial_E .* getfield.(props, :Iz0)
    torsional = mat.G .* (4 .* getfield.(props, :area_mid) .^ 2 .*
                          getfield.(stations, :thickness) ./
                          (2 .* (getfield.(props, :width) .+ getfield.(props, :height))))
    p = Plots.plot(span, bending; label = "EI", xlabel = "Span (m)",
                   ylabel = "Stiffness", title = "Wing spar stiffness")
    Plots.plot!(p, span, torsional; label = "GJ (thin-wall estimate)")
    return p
end

"""
    plot_box_section(prop; title="Box section") -> Plots.Plot

Draw the cross-section outline of a reinforced box section.
"""
function plot_box_section(prop::SectionProperties; title::String = "Box section")
    vertices = vcat(prop.vertices, prop.vertices[1:1, :])
    return Plots.plot(vertices[:, 1], vertices[:, 2]; aspect_ratio = :equal,
                      xlabel = "x (m)", ylabel = "y (m)", title, label = "Section")
end

"""
    plot_spar_topview(stations, wing) -> Plots.Plot

Top-view sketch of the wing planform with the spar overlay.
"""
function plot_spar_topview(stations::Vector{Station}, wing::WingGeometry)
    isempty(stations) && throw(ArgumentError("stations must not be empty"))
    span = getfield.(stations, :span)
    leading_edge = [station.chord / 4 + station.chord * station.limits[1] for station in stations]
    trailing_edge = [station.chord / 4 + station.chord * station.limits[2] for station in stations]
    p = Plots.plot(leading_edge, span; label = "Front web", xlabel = "Chordwise x (m)",
                   ylabel = "Span (m)", title = "Wing spar planform")
    Plots.plot!(p, trailing_edge, span; label = "Rear web")
    Plots.plot!(p, leading_edge, -span; label = "")
    Plots.plot!(p, trailing_edge, -span; label = "")
    Plots.xlims!(p, -wing.root_chord / 2, wing.root_chord)
    return p
end
