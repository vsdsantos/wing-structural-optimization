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

Notes:
- Load vectors in `load` must include both halves, ordered left→right.
- Safety factor 1.5 and uncertainty factor 1.15 are applied internally
  (consistent with original Matlab code).
- Strut angle `wing.constraint_angle` (degrees) is converted to radians.
- The aileron index `.i` in HalfLoads is set to the first station beyond
  the taper break (index 1 when not relevant).
"""
function internal_loads_SS(load::AeroLoad,
                            section_props::Vector{SectionProperties},
                            wing::WingGeometry)
    b     = wing.full_span
    d     = wing.constraint_distance
    theta = deg2rad(wing.constraint_angle)

    # Full-span station vector (add zero-lift tips)
    s   = vcat([-b/2], load.x, [b/2])
    s_0 = s .+ b/2      # shift origin to left tip: 0 … b

    # Centroid x-positions (for torsion shear-centre correction)
    # section_props contains right-half stations only; mirror for left half
    xl_half = [p.centroid[1] for p in section_props]
    xl      = vcat(reverse(xl_half), xl_half)   # left then right
    xl      = vcat([xl[1]], xl, [xl[end]])       # add tip values

    # Aerodynamic loads with safety and uncertainty factors
    sf  = 1.5 * 1.15
    w   = vcat([0.0], load.dLdb  .* sf, [0.0])
    wt  = vcat([0.0], load.dMadb .* sf, [0.0])

    # Strut support indices (nearest panel to (b±d)/2)
    _, i1 = findmin(abs.(s_0 .- (b - d)/2))
    _, i2 = findmin(abs.(s_0 .- (b + d)/2))

    # Cumulative integrals
    int_w   = cumtrapz(s_0, w)
    int_ws  = cumtrapz(s_0, w .* s)
    int_wt1 = cumtrapz(s_0[1:i1], wt[1:i1] .+ xl[1:i1] .* w[1:i1])
    int_wt2 = cumtrapz(s_0[i1:i2], wt[i1:i2] .+ xl[i1:i2] .* w[i1:i2])
    int_wt3 = cumtrapz(s_0[i2:end], wt[i2:end] .+ xl[i2:end] .* w[i2:end])

    # Reaction forces (linear system  A·[Ra; Rb] = B)
    A_mat = [sin(theta)             sin(theta)            ;
             cos(theta)*(b-d)/2     cos(theta)*(b+d)/2    ]
    B_vec = [int_w[end], int_ws[end]]
    Ra, Rb = A_mat \ B_vec

    # --- Section 1: left tip → left support ---
    V1  = int_w[1:i1]
    M1  = s_0[1:i1]  .* V1 .- int_ws[1:i1]
    Mt1 = int_wt1
    F1  = zeros(i1)

    # --- Section 2: between supports ---
    V2  = int_w[i1:i2] .- Ra
    M2  = s_0[i1:i2] .* (V2 .+ Ra) .- int_ws[i1:i2] .- Ra .* (s_0[i1:i2] .- (b-d)/2)
    # Torsion moment distributed between the two supports
    if abs(int_wt2[end]) > eps()
        x_tq = trapz(s_0[i1:i2], wt[i1:i2] .* s_0[i1:i2]) /
               trapz(s_0[i1:i2], wt[i1:i2]  .+ eps())
    else
        x_tq = (s_0[i1] + s_0[i2]) / 2
    end
    Mt2 = -int_wt2 .+ int_wt2[end] * ((b+d)/2 - x_tq) / (d + eps())
    F2  = collect(LinRange(Ra*cos(theta), Rb*cos(theta), i2-i1+1))

    # --- Section 3: right support → right tip ---
    V3  = int_w[i2:end] .- Ra .- Rb
    M3  = s_0[i2:end] .* (V3 .+ Ra .+ Rb) .- int_ws[i2:end] .-
          Ra .* (s_0[i2:end] .- (b-d)/2) .- Rb .* (s_0[i2:end] .- (b+d)/2)
    Mt3 = int_wt3 .- int_wt3[end]
    F3  = zeros(length(s_0) - i2 + 1)

    # --- Assemble full-span arrays (handle overlapping boundary points) ---
    V  = vcat(V1[1:i1-1],  [abs(V1[i1] - V2[1])],   V2[2:end-1],
              [abs(V2[end] - V3[1])],  V3[2:end])
    M  = vcat(M1[1:i1-1],  [max(M1[i1], M2[1])],     M2[2:end-1],
              [max(M2[end], M3[1])],   M3[2:end])
    Mt = vcat(Mt1[1:end-1],
              [abs(Mt1[end]) > eps() ?
               abs(Mt1[end])*max(abs(Mt1[end]),abs(Mt2[1]))/Mt1[end] : 0.0],
              Mt2[2:end-1],
              [abs(Mt3[1]) > eps() ?
               abs(Mt3[1]) *max(abs(Mt2[end]),abs(Mt3[1]))/Mt3[1]   : 0.0],
              Mt3[2:end])
    F  = vcat(F1[1:i1-1], F2[i1:i2], F3[i2+1:end])

    # Remove tip placeholders
    s_inner  = s[2:end-1]
    V_inner  = V[2:end-1]
    M_inner  = M[2:end-1]
    Mt_inner = Mt[2:end-1]
    F_inner  = F[2:end-1]

    n = length(s_inner)
    half = n ÷ 2

    # Right half (second half of the symmetric arrays)
    right = HalfLoads(s_inner[half+1:end],  V_inner[half+1:end],
                      M_inner[half+1:end], Mt_inner[half+1:end],
                      F_inner[half+1:end], 1)

    # Left half (first half, spanwise from root outward after flip)
    left = HalfLoads(reverse(-s_inner[1:half]),  reverse(-V_inner[1:half]),
                     reverse( M_inner[1:half]),  reverse(-Mt_inner[1:half]),
                     reverse( F_inner[1:half]),  1)

    return InternalLoads(left, right)
end

# ---------------------------------------------------------------------------
# Fixed-support (FX) stub
# ---------------------------------------------------------------------------

"""
    internal_loads_FX(load, section_props, wing) -> InternalLoads

**Not fully implemented.** The original `internal_loads_FX.m` contains syntax
errors and was a work-in-progress biplane model. This stub raises an error to
prevent silent incorrect results.
"""
function internal_loads_FX(::AeroLoad, ::Vector{SectionProperties}, ::WingGeometry)
    error("internal_loads_FX is not yet implemented. " *
          "Use constraint_type = \"ss\" for simply-supported analysis.")
end

# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

"""
    compute_internal_loads(load, section_props, wing, constraint_type) -> InternalLoads
"""
function compute_internal_loads(load::AeroLoad,
                                 section_props::Vector{SectionProperties},
                                 wing::WingGeometry,
                                 constraint_type::String)
    if constraint_type == "ss"
        return internal_loads_SS(load, section_props, wing)
    elseif constraint_type == "fx"
        return internal_loads_FX(load, section_props, wing)
    else
        error("Unknown constraint_type: \"$constraint_type\"")
    end
end
