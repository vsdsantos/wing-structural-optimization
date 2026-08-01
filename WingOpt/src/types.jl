"""
Core data structures for the wing structural optimization package.
Replaces all Matlab struct definitions; eliminates file-based data passing.
"""

# ---------------------------------------------------------------------------
# Material
# ---------------------------------------------------------------------------

"""Material mechanical properties."""
struct MaterialProperties
    axial_trac::Float64    # axial tensile strength (Pa)
    axial_comp::Float64    # axial compressive strength (Pa)
    axial_E::Float64       # axial Young's modulus (Pa)
    axial_poison::Float64  # axial Poisson's ratio
    trans_trac::Float64    # transverse tensile strength (Pa)
    trans_comp::Float64    # transverse compressive strength (Pa)
    trans_E::Float64       # transverse Young's modulus (Pa)
    trans_poison::Float64  # transverse Poisson's ratio
    G::Float64             # shear modulus (Pa)
    shear::Float64         # shear strength (Pa)
    dens::Float64          # density (kg/m³)
end

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

"""Wing planform and constraint geometry."""
struct WingGeometry
    root_chord::Float64               # root chord (m)
    full_span::Float64                # full span (m)
    taper_ratio::Vector{Float64}      # chord/root_chord at each transition
    taper_position::Float64           # fraction of semi-span at taper break
    sweep::Vector{Float64}            # 1/4-chord sweep (deg) [NOT IMPLEMENTED]
    dihedral::Vector{Float64}         # dihedral angle (deg)  [NOT IMPLEMENTED]
    torsion::Vector{Float64}          # geometric torsion (deg)[NOT IMPLEMENTED]
    foil::Vector{String}              # airfoil names per segment
    constraint_distance::Float64      # strut attachment span (m)
    constraint_angle::Float64         # strut angle (deg)
    chord_transitions::Vector{Float64}# chord at each span break (m)
    span_transitions::Vector{Float64} # span-wise break positions (m, from root)
end

"""Optimization run configuration."""
Base.@kwdef mutable struct RunConfig
    section_type::String = "reinforced_box"
    constraint_type::String = "ss"
    constraint_gap::Float64 = 0.1
    constraint_angle::Float64 = 0.0
    foil_name::String = "s1223"
end

# ---------------------------------------------------------------------------
# Airfoil geometry
# ---------------------------------------------------------------------------

"""
Parsed airfoil data (replaces the .mat cache written by parse_foil.m).
`ex_lin` and `in_lin` are Interpolations.LinearInterpolation objects mapping
fractional chord x ∈ [0,1] to the y-coordinate on the upper/lower surface.
"""
struct FoilData{IU,IL}
    ex_lin::IU              # upper-surface interpolant  y_upper(x)
    in_lin::IL              # lower-surface interpolant  y_lower(x)
    ext::Matrix{Float64}    # upper surface raw points (x increasing)
    int::Matrix{Float64}    # lower surface raw points (x increasing)
    esp::Float64            # evaluation step (10 % of min data spacing)
end

"""Evaluated airfoil geometry at a specific chord/limits section."""
struct FoilGeometry
    upper::Matrix{Float64}  # upper surface points [x y] (N×2)
    lower::Matrix{Float64}  # lower surface points [x y] (N×2)
    extra_length::Float64   # upper arc length
    intra_length::Float64   # lower arc length
    front_length::Float64   # front-web height (gap at a-limit)
    back_length::Float64    # back-web height  (gap at b-limit)
    perimeter::Float64
    area::Float64           # enclosed area (trapz)
    x::Vector{Float64}      # all x coords (for visualisation)
    y::Vector{Float64}      # all y coords
end

# ---------------------------------------------------------------------------
# Wing stations
# ---------------------------------------------------------------------------

"""One spanwise cross-section station."""
struct Station
    chord::Float64              # local chord (m)
    thickness::Float64          # box wall thickness (m)
    reinforcer::Float64         # reinforcer thickness (m)
    span::Float64               # spanwise position (m from root)
    limits::Vector{Float64}     # [front, back] in fractional chord [0, 1]
end

# ---------------------------------------------------------------------------
# Cross-section structural properties
# ---------------------------------------------------------------------------

"""Structural cross-section properties (replaces all Matlab prop structs)."""
struct SectionProperties
    I::Matrix{Float64}           # inertia tensor (2×2)
    I_max::Matrix{Float64}       # principal inertia tensor (2×2)
    Iz0::Float64                 # I_max[1,1] – used for bending
    Iy0::Float64                 # I_max[2,2]
    theta_max::Float64           # principal axis angle (rad)
    centroid::Vector{Float64}    # centroid [x, y] (m)
    area_mid::Float64            # mid-line enclosed area (m²)
    area::Float64                # material cross-sectional area (m²)
    width::Float64               # section width  (m)
    height::Float64              # section height (m)
    web_length::Float64          # clear web length (m)
    vertices::Matrix{Float64}    # outer corners A,B,C,D (4×2)
    vertices_mid::Matrix{Float64}# mid-line corners (4×2)
    mc_centroid::Vector{Float64} # compression-flange centroid
    mt_centroid::Vector{Float64} # tension-flange centroid
end

# ---------------------------------------------------------------------------
# Parametric geometry
# ---------------------------------------------------------------------------

"""Spar geometry parametrisation over the semi-span."""
struct StructuredGeometry
    x_position::Vector{Float64} # front-web x (m, offset from c/4); one per span break
    width::Vector{Float64}       # spar width (m); one per span break
    thickness::Vector{Float64}   # box wall thickness (m)
    reinforcer::Vector{Float64}  # reinforcer thickness (m)
end

# ---------------------------------------------------------------------------
# Loads
# ---------------------------------------------------------------------------

"""Aerodynamic load distribution (one flight case)."""
struct AeroLoad
    x::Vector{Float64}       # panel centre spanwise positions (m), full span
    dLdb::Vector{Float64}    # lift per unit span (N/m), full span
    dMadb::Vector{Float64}   # pitching moment per unit span (N), full span
end

"""Internal structural loads on one wing half."""
struct HalfLoads
    s::Vector{Float64}   # spanwise stations (m)
    V::Vector{Float64}   # shear force (N)
    M::Vector{Float64}   # bending moment (N·m)
    Mt::Vector{Float64}  # torsional moment (N·m)
    F::Vector{Float64}   # axial force (N)
    i::Int               # index at which aileron region starts (for Report-45)
end

"""Internal loads for the full wing (both halves)."""
struct InternalLoads
    L::HalfLoads
    R::HalfLoads
end

# ---------------------------------------------------------------------------
# Optimisation
# ---------------------------------------------------------------------------

"""Result returned by every optimizer wrapper."""
struct OptResult
    x::Vector{Float64}
    fval::Float64
    exitflag::Int
    output::Any
end

"""
All data produced by a single fitness-function evaluation.
Replaces the ad-hoc `data` struct assembled at the end of fit_function.m.
"""
mutable struct FitData
    stations::Vector{Station}
    section_props::Vector{SectionProperties}
    wing::WingGeometry
    mass::Float64
    dm::Vector{Float64}
    min_ms::Float64
    ms::Vector{Float64}
    internal_loads::Vector{InternalLoads}
end
