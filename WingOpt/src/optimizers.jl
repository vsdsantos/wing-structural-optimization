"""
Optimiser wrappers.
Translates `optimizers/`: ga_opt.m, ps_opt.m, psw_opt.m, sa_opt.m.

Uses Metaheuristics.jl for GA / PSO / SA and Optim.jl for pattern search.
"""

# ---------------------------------------------------------------------------
# Shared bounds and constraint definitions
# ---------------------------------------------------------------------------

"""Default optimisation bounds matching `start_opt.m`."""
const DEFAULT_LB = [0.01, 0.01, 0.01, 1.0, 0.002]
const DEFAULT_UB = [0.60, 0.40, 0.40, 4.0, 0.010]

"""Linear inequality matrix from `start_opt.m`: [1 1 0 0 0; 1 0 1 0 0]·x ≤ [0.75; 0.75]."""
const AINEQ = [
    1.0 1.0 0.0 0.0 0.0;
    1.0 0.0 1.0 0.0 0.0
]
const BINEQ = [0.75, 0.75]

# ---------------------------------------------------------------------------
# Genetic algorithm  (ga_opt.m)
# ---------------------------------------------------------------------------

"""
    run_ga(fitness_fn; lb, ub, max_generations, population_size,
           parallel, seed) -> OptResult

Genetic algorithm optimiser using Metaheuristics.GA.
"""
function run_ga(
    fitness_fn::Function;
    lb::Vector{Float64} = DEFAULT_LB,
    ub::Vector{Float64} = DEFAULT_UB,
    max_generations::Int = 100,
    population_size::Int = 200,
    parallel::Bool = false,
    seed::Int = 42,
)
    error("Not implemented")
end

# ---------------------------------------------------------------------------
# Particle swarm  (psw_opt.m)
# ---------------------------------------------------------------------------

"""
    run_pso(fitness_fn; lb, ub, max_iterations, swarm_size,
            parallel, seed) -> OptResult

Particle swarm optimiser using Metaheuristics.PSO.
"""
function run_pso(
    fitness_fn::Function;
    lb::Vector{Float64} = DEFAULT_LB,
    ub::Vector{Float64} = DEFAULT_UB,
    max_iterations::Int = 200,
    swarm_size::Int = 100,
    parallel::Bool = false,
    seed::Int = 42,
)
    error("Not implemented")
end

# ---------------------------------------------------------------------------
# Simulated annealing  (sa_opt.m)
# ---------------------------------------------------------------------------

"""
    run_sa(fitness_fn, x0; lb, ub, max_iterations, seed) -> OptResult

Simulated annealing using Metaheuristics.SA.
"""
function run_sa(
    fitness_fn::Function,
    x0::Vector{Float64};
    lb::Vector{Float64} = DEFAULT_LB,
    ub::Vector{Float64} = DEFAULT_UB,
    max_iterations::Int = 5000,
    seed::Int = 42,
)
    error("Not implemented")
end

# ---------------------------------------------------------------------------
# Pattern search / NelderMead  (ps_opt.m)
# ---------------------------------------------------------------------------

"""
    run_ps(fitness_fn, x0; lb, ub, max_iterations) -> OptResult

Pattern search (Nelder-Mead) using Optim.jl.
"""
function run_ps(
    fitness_fn::Function,
    x0::Vector{Float64};
    lb::Vector{Float64} = DEFAULT_LB,
    ub::Vector{Float64} = DEFAULT_UB,
    max_iterations::Int = 1000,
)
    error("Not implemented")
end

# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

"""
    run_optimizer(optimizer, fitness_fn; x0, lb, ub, kwargs...) -> OptResult

Dispatch to the named optimiser.

`optimizer` – one of `:ga`, `:pso`, `:sa`, `:ps`
`x0`        – initial point (required for `:sa` and `:ps`)
"""
function run_optimizer(
    optimizer::Symbol,
    fitness_fn::Function;
    x0::Union{Nothing,Vector{Float64}} = nothing,
    lb::Vector{Float64} = DEFAULT_LB,
    ub::Vector{Float64} = DEFAULT_UB,
    kwargs...,
)
    error("Not implemented")
end
