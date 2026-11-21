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
    run_case == 1 ||
        throw(ArgumentError("unsupported run case $run_case; only run case 1 is defined"))
    return RunConfig()
end

"""
    build_wing_geometry(cfg::RunConfig) -> WingGeometry

Build the wing geometry struct from the run configuration.
Mirrors the values in `wing_geometry.m`.
"""
function build_wing_geometry(cfg::RunConfig)
    cfg.constraint_gap >= 0 || throw(ArgumentError("constraint_gap must be nonnegative"))
    return WingGeometry(
        1.0,
        7.0,
        [1.0, 0.5],
        0.75,
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0],
        [cfg.foil_name, cfg.foil_name],
        cfg.constraint_gap,
        cfg.constraint_angle,
        [1.0, 1.0, 0.5],
        [0.0, 2.625, 3.5],
    )
end
