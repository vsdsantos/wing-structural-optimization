"""
Main entry point.
Translates `start_opt.m` from a workspace-mutation script to a clean function.
"""

"""
    run_optimization(run_case=1;
                     optimizer=:ga,
                     parallel=false,
                     foil_dir="foils",
                     loads=nothing,
                     n_panels=20,
                     kwargs...) -> (OptResult, FitData)

Full optimization pipeline:

1. Build run configuration and wing geometry.
2. Parse airfoil data.
3. Build (or use supplied) aerodynamic loads.
4. Run the selected optimiser.
5. Evaluate the best design and return `(OptResult, FitData)`.

Arguments:
- `run_case`  – integer run-case selector (see `config.jl`)
- `optimizer` – `:ga` | `:pso` | `:sa` | `:ps`
- `parallel`  – enable parallel function evaluations (currently a hint to the solver)
- `foil_dir`  – directory containing `.dat` airfoil files
- `loads`     – `Vector{AeroLoad}` override; if `nothing`, an elliptic example is used
- `n_panels`  – number of aerodynamic panels when using the example load
- `kwargs`    – forwarded to the chosen optimiser (e.g. `max_generations=50`)

Optimisation variables (1-based):
    X[1] – web_pos_tip       ∈ [0.01, 0.60]  (fraction of chord)
    X[2] – width_tip         ∈ [0.01, 0.40]
    X[3] – width_root        ∈ [0.01, 0.40]
    X[4] – box_thickness     ∈ {1,2,3,4}     (integer, 1 mm units)
    X[5] – reinforcer_thick  ∈ [0.002, 0.010] (m)

Constraints:
    X[1] + X[2] ≤ 0.75   (tip spar must fit within chord)
    X[1] + X[3] ≤ 0.75   (root spar must fit within chord)
"""
function run_optimization(run_case::Int=1;
                           optimizer::Symbol=:ga,
                           parallel::Bool=false,
                           foil_dir::String="foils",
                           loads::Union{Nothing,Vector{AeroLoad}}=nothing,
                           n_panels::Int=20,
                           kwargs...)
    error("Not implemented")
end
