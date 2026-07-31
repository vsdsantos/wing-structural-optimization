"""
Optimiser wrappers.
Translates `optimizers/`: ga_opt.m, ps_opt.m, psw_opt.m, sa_opt.m.

Uses Metaheuristics.jl for GA / PSO / SA and Optim.jl for pattern search.

Linear inequality constraints `Aineq·x ≤ bineq` are enforced via a penalty
added to the objective (both Metaheuristics and Optim have limited native
support for linear constraints in the general black-box setting).
"""

# ---------------------------------------------------------------------------
# Shared bounds and constraint definitions
# ---------------------------------------------------------------------------

"""Default optimisation bounds matching `start_opt.m`."""
const DEFAULT_LB = [0.01, 0.01, 0.01, 1.0,  0.002]
const DEFAULT_UB = [0.60, 0.40, 0.40, 4.0,  0.010]

"""Linear inequality matrix from `start_opt.m`: [1 1 0 0 0; 1 0 1 0 0]·x ≤ [0.75; 0.75]."""
const AINEQ = [1.0 1.0 0.0 0.0 0.0;
               1.0 0.0 1.0 0.0 0.0]
const BINEQ = [0.75, 0.75]

"""
    penalty(x, Aineq, bineq; mu=1e6) -> Float64

Return a penalty value for linear constraint violations `Aineq·x ≤ bineq`.
"""
function penalty(x::AbstractVector,
                  Aineq::Matrix{Float64}=AINEQ,
                  bineq::Vector{Float64}=BINEQ;
                  mu::Float64=1e6)
    violation = max.(Aineq * x .- bineq, 0.0)
    return mu * sum(violation.^2)
end

# ---------------------------------------------------------------------------
# Objective wrapper (with bounds projection and constraint penalty)
# ---------------------------------------------------------------------------

"""
    make_objective(fitness_fn, lb, ub) -> Function

Wrap the fitness function to:
1. Clamp `x` within bounds.
2. Round `x[4]` to the nearest integer (integer variable).
3. Add a penalty for linear constraint violations.
"""
function make_objective(fitness_fn::Function,
                         lb::Vector{Float64}=DEFAULT_LB,
                         ub::Vector{Float64}=DEFAULT_UB)
    function obj(x)
        x_c = clamp.(x, lb, ub)
        x_c[4] = round(x_c[4])
        raw, _ = fitness_fn(x_c)
        return raw + penalty(x_c)
    end
    return obj
end

# ---------------------------------------------------------------------------
# Genetic algorithm  (ga_opt.m)
# ---------------------------------------------------------------------------

"""
    run_ga(fitness_fn; lb, ub, max_generations, population_size,
           parallel, seed) -> OptResult

Genetic algorithm optimiser using Metaheuristics.GA.
"""
function run_ga(fitness_fn::Function;
                lb::Vector{Float64}=DEFAULT_LB,
                ub::Vector{Float64}=DEFAULT_UB,
                max_generations::Int=100,
                population_size::Int=200,
                parallel::Bool=false,
                seed::Int=42)
    obj = make_objective(fitness_fn, lb, ub)

    bounds = boxconstraints(lb=lb, ub=ub)
    options = Options(seed=seed,
                      iterations=max_generations,
                      f_calls_limit=max_generations * population_size * 10)

    alg = GA(;
        N         = population_size,
        options   = options
    )

    result = optimize(obj, bounds, alg)

    x_best = minimizer(result)
    x_best[4] = round(x_best[4])
    return OptResult(x_best, minimum(result), 0, result)
end

# ---------------------------------------------------------------------------
# Particle swarm  (psw_opt.m)
# ---------------------------------------------------------------------------

"""
    run_pso(fitness_fn; lb, ub, max_iterations, swarm_size,
            parallel, seed) -> OptResult

Particle swarm optimiser using Metaheuristics.PSO.
"""
function run_pso(fitness_fn::Function;
                 lb::Vector{Float64}=DEFAULT_LB,
                 ub::Vector{Float64}=DEFAULT_UB,
                 max_iterations::Int=200,
                 swarm_size::Int=100,
                 parallel::Bool=false,
                 seed::Int=42)
    obj = make_objective(fitness_fn, lb, ub)
    bounds = boxconstraints(lb=lb, ub=ub)
    options = Options(seed=seed, iterations=max_iterations)
    alg = PSO(; N=swarm_size, options=options)

    result = optimize(obj, bounds, alg)
    x_best = minimizer(result);  x_best[4] = round(x_best[4])
    return OptResult(x_best, minimum(result), 0, result)
end

# ---------------------------------------------------------------------------
# Simulated annealing  (sa_opt.m)
# ---------------------------------------------------------------------------

"""
    run_sa(fitness_fn, x0; lb, ub, max_iterations, seed) -> OptResult

Simulated annealing using Metaheuristics.SA.
"""
function run_sa(fitness_fn::Function,
                x0::Vector{Float64};
                lb::Vector{Float64}=DEFAULT_LB,
                ub::Vector{Float64}=DEFAULT_UB,
                max_iterations::Int=5000,
                seed::Int=42)
    obj = make_objective(fitness_fn, lb, ub)
    bounds = boxconstraints(lb=lb, ub=ub)
    options = Options(seed=seed, iterations=max_iterations)
    alg = SA(; options=options)

    result = optimize(obj, bounds, alg)
    x_best = minimizer(result);  x_best[4] = round(x_best[4])
    return OptResult(x_best, minimum(result), 0, result)
end

# ---------------------------------------------------------------------------
# Pattern search / NelderMead  (ps_opt.m)
# ---------------------------------------------------------------------------

"""
    run_ps(fitness_fn, x0; lb, ub, max_iterations) -> OptResult

Pattern search (Nelder-Mead) using Optim.jl.
"""
function run_ps(fitness_fn::Function,
                x0::Vector{Float64};
                lb::Vector{Float64}=DEFAULT_LB,
                ub::Vector{Float64}=DEFAULT_UB,
                max_iterations::Int=1000)
    obj = make_objective(fitness_fn, lb, ub)

    result = Optim.optimize(obj, lb, ub, x0,
                            Optim.Fminbox(Optim.NelderMead()),
                            Optim.Options(iterations=max_iterations))
    x_best = Optim.minimizer(result);  x_best[4] = round(x_best[4])
    return OptResult(x_best, Optim.minimum(result),
                     Optim.converged(result) ? 1 : 0, result)
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
function run_optimizer(optimizer::Symbol,
                        fitness_fn::Function;
                        x0::Union{Nothing,Vector{Float64}}=nothing,
                        lb::Vector{Float64}=DEFAULT_LB,
                        ub::Vector{Float64}=DEFAULT_UB,
                        kwargs...)
    if optimizer == :ga
        return run_ga(fitness_fn; lb=lb, ub=ub, kwargs...)
    elseif optimizer == :pso
        return run_pso(fitness_fn; lb=lb, ub=ub, kwargs...)
    elseif optimizer == :sa
        x0 === nothing && error("run_optimizer: x0 required for :sa")
        return run_sa(fitness_fn, x0; lb=lb, ub=ub, kwargs...)
    elseif optimizer == :ps
        x0 === nothing && error("run_optimizer: x0 required for :ps")
        return run_ps(fitness_fn, x0; lb=lb, ub=ub, kwargs...)
    else
        error("Unknown optimizer: $optimizer. Use :ga, :pso, :sa, or :ps.")
    end
end
