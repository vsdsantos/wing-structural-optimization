"""
Configuration builders.
Translates `run_cases.m` + `wing_geometry.m` from script-based
workspace mutation to clean, testable functions.
"""

"""
    build_run_config(run_case::Int) -> RunConfig

Return the RunConfig for the requested run case number.
Add new cases by extending the if/elseif chain below.
"""
function build_run_config(run_case::Int)
    if run_case == 1
        return RunConfig(
            section_type    = "reinforced_box",
            constraint_type = "ss",
            constraint_gap  = 0.1,
            constraint_angle = 0.0,
            foil_name       = "s1223"
        )
    else
        error("Unknown run case: $run_case. Define it in config.jl.")
    end
end

"""
    build_wing_geometry(cfg::RunConfig) -> WingGeometry

Build the wing geometry struct from the run configuration.
Mirrors the values in `wing_geometry.m`.

The planform is a double-tapered wing:

```
       root_chord
    |--------------|
    |              |---
    |              |  | pci·b/2
    |              |  |
     \\            /
      \\          /
       \\________/
```
"""
function build_wing_geometry(cfg::RunConfig)
    root_chord      = 1.0          # m
    full_span       = 7.0          # m
    taper_ratio     = [1.0, 0.5]   # chord / root_chord
    taper_position  = 0.75         # fraction of semi-span at taper break

    chord_transitions = root_chord .* vcat([1.0], taper_ratio)   # [1.0, 1.0, 0.5]
    half_span = full_span / 2
    span_transitions  = vcat([0.0],
                              cumsum(half_span .* [taper_position,
                                                   1.0 - taper_position]))

    return WingGeometry(
        root_chord,
        full_span,
        taper_ratio,
        taper_position,
        [0.0, 0.0],      # sweep   (NOT IMPLEMENTED)
        [0.0, 0.0],      # dihedral(NOT IMPLEMENTED)
        [0.0, 0.0],      # torsion (NOT IMPLEMENTED)
        [cfg.foil_name, cfg.foil_name],
        cfg.constraint_gap,
        cfg.constraint_angle,
        chord_transitions,
        span_transitions
    )
end
