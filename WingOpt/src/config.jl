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
    error("Not implemented")
end

"""
    build_wing_geometry(cfg::RunConfig) -> WingGeometry

Build the wing geometry struct from the run configuration.
Mirrors the values in `wing_geometry.m`.
"""
function build_wing_geometry(cfg::RunConfig)
    error("Not implemented")
end
