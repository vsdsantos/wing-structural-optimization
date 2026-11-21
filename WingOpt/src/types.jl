"""
Core data structures for the wing structural optimization package.
Replaces all Matlab struct definitions; eliminates file-based data passing.
"""

# ---------------------------------------------------------------------------
# Material
# ---------------------------------------------------------------------------

"""
Material mechanical properties in SI units: strengths/moduli in Pa, density in
kg/m³, and Poisson ratios dimensionless.
"""
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
    function MaterialProperties(values::Vararg{Real,11})
        v = Float64.(values)
        all(isfinite, v) || throw(ArgumentError("material properties must be finite"))
        all(>(0), v[[1, 2, 3, 5, 6, 7, 9, 10, 11]]) ||
            throw(DomainError(v, "strengths, moduli, and density must be positive"))
        all(ν -> 0 <= ν < 0.5, v[[4, 8]]) ||
            throw(DomainError(v[[4, 8]], "Poisson ratios must be in [0, 0.5)"))
        return new(v...)
    end
end

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

"""
Wing planform and constraint geometry. Lengths are in metres, sweep/dihedral/
torsion and constraint angles are in degrees, and span coordinates increase
outward from the root. All geometric values are SI-based.
"""
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
    function WingGeometry(
        root_chord::Real,
        full_span::Real,
        taper_ratio::AbstractVector{<:Real},
        taper_position::Real,
        sweep::AbstractVector{<:Real},
        dihedral::AbstractVector{<:Real},
        torsion::AbstractVector{<:Real},
        foil::AbstractVector{<:AbstractString},
        constraint_distance::Real,
        constraint_angle::Real,
        chord_transitions::AbstractVector{<:Real},
        span_transitions::AbstractVector{<:Real},
    )
        scalar_values =
            (root_chord, full_span, taper_position, constraint_distance, constraint_angle)
        vectors =
            (taper_ratio, sweep, dihedral, torsion, chord_transitions, span_transitions)
        all(isfinite, scalar_values) && all(v -> all(isfinite, v), vectors) ||
            throw(ArgumentError("wing geometry values must be finite"))
        root_chord > 0 && full_span > 0 || throw(
            DomainError(
                (root_chord, full_span),
                "root chord and full span must be positive",
            ),
        )
        0 < taper_position <= 1 ||
            throw(DomainError(taper_position, "taper position must be in (0, 1]"))
        all(>(0), taper_ratio) ||
            throw(DomainError(taper_ratio, "taper ratios must be positive"))
        n = length(span_transitions)
        n >= 2 && length(chord_transitions) == n || throw(
            DimensionMismatch("at least two matching chord/span transitions are required"),
        )
        length(taper_ratio) == n - 1 &&
        length(sweep) == n - 1 &&
        length(dihedral) == n - 1 &&
        length(torsion) == n - 1 &&
        length(foil) == n - 1 || throw(
            DimensionMismatch(
                "segment data must contain one entry per transition interval",
            ),
        )
        all(>(0), chord_transitions) ||
            throw(DomainError(chord_transitions, "chord lengths must be positive"))
        all(diff(span_transitions) .> 0) ||
            throw(ArgumentError("span transitions must be strictly increasing"))
        isapprox(first(span_transitions), 0; atol = 1e-10) &&
        isapprox(last(span_transitions), full_span / 2; rtol = 1e-10, atol = 1e-10) ||
            throw(
                DomainError(
                    span_transitions,
                    "span transitions must run from root to semi-span tip",
                ),
            )
        0 <= constraint_distance <= full_span / 2 || throw(
            DomainError(
                constraint_distance,
                "constraint distance must lie on the semi-span",
            ),
        )
        all(name -> !isempty(strip(name)), foil) ||
            throw(ArgumentError("airfoil names must not be empty"))
        return new(
            Float64(root_chord),
            Float64(full_span),
            Float64.(taper_ratio),
            Float64(taper_position),
            Float64.(sweep),
            Float64.(dihedral),
            Float64.(torsion),
            String.(foil),
            Float64(constraint_distance),
            Float64(constraint_angle),
            Float64.(chord_transitions),
            Float64.(span_transitions),
        )
    end
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
    function FoilData(
        ex_lin::IU,
        in_lin::IL,
        ext::AbstractMatrix{<:Real},
        int::AbstractMatrix{<:Real},
        esp::Real,
    ) where {IU,IL}
        for (surface, label) in ((ext, "upper"), (int, "lower"))
            size(surface, 2) == 2 && size(surface, 1) >= 2 || throw(
                DimensionMismatch(
                    "$label airfoil surface must contain at least two [x, y] points",
                ),
            )
            all(isfinite, surface) ||
                throw(ArgumentError("$label airfoil coordinates must be finite"))
            all(diff(surface[:, 1]) .> 0) || throw(
                ArgumentError("$label airfoil x-coordinates must be strictly increasing"),
            )
            0 <= first(surface[:, 1]) && last(surface[:, 1]) <= 1 || throw(
                DomainError(surface[:, 1], "airfoil x-coordinates must lie in [0, 1]"),
            )
        end
        max(first(ext[:, 1]), first(int[:, 1])) < min(last(ext[:, 1]), last(int[:, 1])) ||
            throw(
                ArgumentError(
                    "upper and lower airfoil surfaces must have overlapping x-ranges",
                ),
            )
        isfinite(esp) && esp > 0 || throw(
            DomainError(esp, "airfoil evaluation spacing must be positive and finite"),
        )
        return new{IU,IL}(
            ex_lin,
            in_lin,
            Matrix{Float64}(ext),
            Matrix{Float64}(int),
            Float64(esp),
        )
    end
end

"""Evaluated airfoil geometry at a specific section; all coordinates and lengths are metres."""
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
    function FoilGeometry(
        upper::AbstractMatrix{<:Real},
        lower::AbstractMatrix{<:Real},
        extra_length::Real,
        intra_length::Real,
        front_length::Real,
        back_length::Real,
        perimeter::Real,
        area::Real,
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
    )
        for (surface, label) in ((upper, "upper"), (lower, "lower"))
            size(surface, 2) == 2 && size(surface, 1) >= 2 || throw(
                DimensionMismatch(
                    "$label foil surface must have at least two [x, y] points",
                ),
            )
            all(isfinite, surface) ||
                throw(ArgumentError("$label foil surface coordinates must be finite"))
            all(diff(surface[:, 1]) .> 0) ||
                throw(ArgumentError("$label foil surface x-coordinates must increase"))
        end
        length(x) == length(y) ||
            throw(DimensionMismatch("foil display x and y arrays must match"))
        all(isfinite, x) && all(isfinite, y) ||
            throw(ArgumentError("foil display coordinates must be finite"))
        metrics = (extra_length, intra_length, front_length, back_length, perimeter, area)
        all(isfinite, metrics) ||
            throw(ArgumentError("foil geometry metrics must be finite"))
        extra_length > 0 &&
        intra_length > 0 &&
        front_length >= 0 &&
        back_length >= 0 &&
        perimeter > 0 &&
        area > 0 || throw(DomainError(metrics, "foil lengths and area must be physical"))
        return new(
            Matrix{Float64}(upper),
            Matrix{Float64}(lower),
            Float64(extra_length),
            Float64(intra_length),
            Float64(front_length),
            Float64(back_length),
            Float64(perimeter),
            Float64(area),
            Float64.(x),
            Float64.(y),
        )
    end
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
    function Station(
        chord::Real,
        thickness::Real,
        reinforcer::Real,
        span::Real,
        limits::AbstractVector{<:Real},
    )
        all(isfinite, (chord, thickness, reinforcer, span)) && all(isfinite, limits) ||
            throw(ArgumentError("station values must be finite"))
        chord > 0 && thickness > 0 && reinforcer >= 0 && span >= 0 || throw(
            DomainError(
                (chord, thickness, reinforcer, span),
                "station dimensions are nonphysical",
            ),
        )
        length(limits) == 2 ||
            throw(DimensionMismatch("station limits must be [front, back]"))
        0 <= limits[1] < limits[2] <= 1 ||
            throw(DomainError(limits, "station limits must satisfy 0 ≤ front < back ≤ 1"))
        return new(
            Float64(chord),
            Float64(thickness),
            Float64(reinforcer),
            Float64(span),
            Float64.(limits),
        )
    end
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
    function SectionProperties(
        I::AbstractMatrix{<:Real},
        I_max::AbstractMatrix{<:Real},
        Iz0::Real,
        Iy0::Real,
        theta_max::Real,
        centroid::AbstractVector{<:Real},
        area_mid::Real,
        area::Real,
        width::Real,
        height::Real,
        web_length::Real,
        vertices::AbstractMatrix{<:Real},
        vertices_mid::AbstractMatrix{<:Real},
        mc_centroid::AbstractVector{<:Real},
        mt_centroid::AbstractVector{<:Real},
    )
        size(I) == size(I_max) == (2, 2) ||
            throw(DimensionMismatch("section inertia tensors must be 2×2"))
        size(vertices) == size(vertices_mid) == (4, 2) ||
            throw(DimensionMismatch("section vertices must be 4×2"))
        length(centroid) == length(mc_centroid) == length(mt_centroid) == 2 ||
            throw(DimensionMismatch("section centroids must contain [x, y]"))
        scalars = (Iz0, Iy0, theta_max, area_mid, area, width, height, web_length)
        arrays = (I, I_max, centroid, vertices, vertices_mid, mc_centroid, mt_centroid)
        all(isfinite, scalars) && all(v -> all(isfinite, v), arrays) ||
            throw(ArgumentError("section properties must be finite"))
        all(>(0), (Iz0, Iy0, area_mid, area, width, height)) && web_length >= 0 || throw(
            DomainError(
                scalars,
                "section dimensions, areas, and inertias must be physical",
            ),
        )
        isapprox(I, I'; rtol = 1e-10, atol = 1e-14) &&
        isapprox(I_max, I_max'; rtol = 1e-10, atol = 1e-14) ||
            throw(ArgumentError("section inertia tensors must be symmetric"))
        return new(
            Matrix{Float64}(I),
            Matrix{Float64}(I_max),
            Float64(Iz0),
            Float64(Iy0),
            Float64(theta_max),
            Float64.(centroid),
            Float64(area_mid),
            Float64(area),
            Float64(width),
            Float64(height),
            Float64(web_length),
            Matrix{Float64}(vertices),
            Matrix{Float64}(vertices_mid),
            Float64.(mc_centroid),
            Float64.(mt_centroid),
        )
    end
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
    function StructuredGeometry(
        x_position::AbstractVector{<:Real},
        width::AbstractVector{<:Real},
        thickness::AbstractVector{<:Real},
        reinforcer::AbstractVector{<:Real},
    )
        n = length(x_position)
        all(length(v) == n for v in (width, thickness, reinforcer)) || throw(
            DimensionMismatch("all structural geometry arrays must have equal length"),
        )
        all(v -> all(isfinite, v), (x_position, width, thickness, reinforcer)) ||
            throw(ArgumentError("structural geometry values must be finite"))
        all(>(0), width) && all(>(0), thickness) && all(>=(0), reinforcer) || throw(
            DomainError(
                (width, thickness, reinforcer),
                "width/thickness must be positive and reinforcer nonnegative",
            ),
        )
        return new(
            Float64.(x_position),
            Float64.(width),
            Float64.(thickness),
            Float64.(reinforcer),
        )
    end
end

# ---------------------------------------------------------------------------
# Loads
# ---------------------------------------------------------------------------

"""
Aerodynamic load distribution (one flight case). Coordinates are metres over
the full span; lift is positive upward (N/m), and pitching moment is positive
nose-up (N per unit span).
"""
struct AeroLoad
    x::Vector{Float64}       # panel centre spanwise positions (m), full span
    dLdb::Vector{Float64}    # lift per unit span (N/m), full span
    dMadb::Vector{Float64}   # pitching moment per unit span (N), full span
    function AeroLoad(
        x::AbstractVector{<:Real},
        dLdb::AbstractVector{<:Real},
        dMadb::AbstractVector{<:Real},
    )
        length(x) == length(dLdb) == length(dMadb) || throw(
            DimensionMismatch(
                "load coordinates, lift, and pitching moment must have equal lengths",
            ),
        )
        !isempty(x) ||
            throw(ArgumentError("aerodynamic load must contain at least one sample"))
        all(v -> all(isfinite, v), (x, dLdb, dMadb)) ||
            throw(ArgumentError("aerodynamic load values must be finite"))
        all(diff(x) .> 0) ||
            length(x) == 1 ||
            throw(ArgumentError("load coordinates must be strictly increasing"))
        return new(Float64.(x), Float64.(dLdb), Float64.(dMadb))
    end
end

"""
Internal structural loads on one wing half. `s` is in metres from the root;
forces are N and moments are N·m. `F` is positive in tension, with resultant
signs otherwise taken on the positive-span section face (right-hand rule).
"""
struct HalfLoads
    s::Vector{Float64}   # spanwise stations (m, increasing from root to tip)
    V::Vector{Float64}   # shear force (N)
    M::Vector{Float64}   # bending moment (N·m)
    Mt::Vector{Float64}  # torsional moment (N·m)
    F::Vector{Float64}   # axial force (N, positive in tension)
    i::Int               # index at which aileron region starts (for Report-45)
    function HalfLoads(
        s::AbstractVector{<:Real},
        V::AbstractVector{<:Real},
        M::AbstractVector{<:Real},
        Mt::AbstractVector{<:Real},
        F::AbstractVector{<:Real},
        i::Integer,
    )
        n = length(s)
        all(length(v) == n for v in (V, M, Mt, F)) ||
            throw(DimensionMismatch("all internal load arrays must have equal lengths"))
        all(v -> all(isfinite, v), (s, V, M, Mt, F)) ||
            throw(ArgumentError("internal loads must be finite"))
        all(diff(s) .>= 0) || throw(ArgumentError("span coordinates must be nondecreasing"))
        1 <= i <= n + 1 || throw(BoundsError(s, i))
        return new(Float64.(s), Float64.(V), Float64.(M), Float64.(Mt), Float64.(F), Int(i))
    end
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
