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
"""
function internal_loads_SS(load::AeroLoad,
                            section_props::Vector{SectionProperties},
                            wing::WingGeometry)
    error("Not implemented")
end

# ---------------------------------------------------------------------------
# Fixed-support (FX) stub
# ---------------------------------------------------------------------------

"""
    internal_loads_FX(load, section_props, wing) -> InternalLoads

**Not yet implemented.** The original `internal_loads_FX.m` contains syntax
errors and was a work-in-progress biplane model.
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
    error("Not implemented")
end
