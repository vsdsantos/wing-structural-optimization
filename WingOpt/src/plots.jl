"""
Visualization functions.
Translates `plots/`: plot_esforcos.m, plot_longarina.m,
plot_prop_geo.m, plot_sec_retangular.m (reinforced box), plot_sec_perfilada.m.

All functions return the Plots.jl plot object so the caller can save/display.
"""

# ---------------------------------------------------------------------------
# Internal load diagrams  (plot_esforcos.m)
# ---------------------------------------------------------------------------

"""
    plot_internal_loads(loads_half; title_prefix="Right wing") -> Plots.Plot

Plot shear force, bending moment, torsional moment, and axial force
diagrams for one wing half.
"""
function plot_internal_loads(hl::HalfLoads; title_prefix::String="Right wing")
    p1 = Plots.plot(hl.s, hl.V;  label="", ylabel="V (N)",    title="Shear")
    p2 = Plots.plot(hl.s, hl.M;  label="", ylabel="M (N·m)",  title="Bending moment")
    p3 = Plots.plot(hl.s, hl.Mt; label="", ylabel="Mt (N·m)", title="Torsional moment")
    p4 = Plots.plot(hl.s, hl.F;  label="", ylabel="F (N)",    title="Axial force")

    for p in (p1, p2, p3, p4)
        Plots.xlabel!(p, "span (m)")
        Plots.yaxis!(p, grid=true)
    end

    return Plots.plot(p1, p2, p3, p4; layout=(2,2),
                      plot_title=title_prefix * " – Internal loads")
end

# ---------------------------------------------------------------------------
# Section properties  (plot_prop_geo.m)
# ---------------------------------------------------------------------------

"""
    plot_stiffness(stations, props, mat) -> Plots.Plot

Plot bending stiffness EI and torsional stiffness GJ vs semi-span.
"""
function plot_stiffness(stations::Vector{Station},
                         props::Vector{SectionProperties},
                         mat::MaterialProperties)
    s  = [st.span for st in stations]
    I  = [p.Iz0 for p in props]
    J  = [p.Iz0 + p.Iy0 for p in props]
    EI = mat.axial_E .* I
    GJ = mat.G       .* J

    p = Plots.plot(s, EI; label="EI", ylabel="Stiffness (N·m²)",
                   xlabel="semi-span (m)", title="EI & GJ vs span")
    Plots.plot!(p, s, GJ; label="GJ")
    return p
end

# ---------------------------------------------------------------------------
# Box section cross-section  (plot_sec_retangular.m)
# ---------------------------------------------------------------------------

"""
    plot_box_section(prop; title="Box section") -> Plots.Plot

Draw the cross-section outline of a reinforced box section.
"""
function plot_box_section(prop::SectionProperties; title::String="Box section")
    V = prop.vertices         # outer rectangle
    Vm = prop.vertices_mid    # inner (mid-line) rectangle

    # Close the rectangles for plotting
    close_rect(R) = vcat(R, R[[1], :])
    Vo  = close_rect(V)
    Vmo = close_rect(Vm)

    p = Plots.plot(Vo[:, 1], Vo[:, 2]; aspect_ratio=:equal,
                   label="outer", lw=2, color=:blue)
    Plots.plot!(p, Vmo[:, 1], Vmo[:, 2]; label="mid-line", lw=1,
                color=:blue, linestyle=:dash)
    Plots.scatter!(p, [prop.centroid[1]], [prop.centroid[2]];
                   label="centroid", markersize=5, color=:red)
    Plots.xlabel!(p, "x (m)");  Plots.ylabel!(p, "y (m)")
    Plots.title!(p, title)
    return p
end

# ---------------------------------------------------------------------------
# Spar top view  (plot_longarina.m)
# ---------------------------------------------------------------------------

"""
    plot_spar_topview(stations, wing) -> Plots.Plot

Top-view sketch of the wing planform with the spar overlay.
"""
function plot_spar_topview(stations::Vector{Station}, wing::WingGeometry)
    sp = wing.span_transitions
    ch = wing.chord_transitions

    p = Plots.plot(; aspect_ratio=:equal, yflip=true,
                   title="Spar – top view", xlabel="span (m)", ylabel="chord (m)")

    # Wing outline (LE and TE)
    for (si, ci) in zip(sp, ch)
        Plots.plot!(p, [si, si], [-ci*0.25, ci*0.75]; color=:black, linestyle=:dash, label="")
    end
    Plots.plot!(p, sp, ch .* 0.75;  color=:black, linestyle=:dash, label="")
    Plots.plot!(p, sp, -ch .* 0.25; color=:black, linestyle=:dash, label="")

    # Spar front and back edges
    s_vals = [st.span  for st in stations]
    a_vals = [st.limits[1] * st.chord - st.chord / 4 for st in stations]
    b_vals = [st.limits[2] * st.chord - st.chord / 4 for st in stations]
    Plots.plot!(p, s_vals, a_vals; color=:blue, label="spar front")
    Plots.plot!(p, s_vals, b_vals; color=:blue, label="spar back")
    return p
end
