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
    population_size >= 2 && max_generations >= 1 ||
        throw(ArgumentError("population_size must be ≥2 and generations must be positive"))
    rng = MersenneTwister(seed)
    population = [_project_design(lb .+ rand(rng, length(lb)) .* (ub .- lb), lb, ub)
                  for _ in 1:population_size]
    scores = [_objective(fitness_fn, x) for x in population]
    best_idx = argmin(scores)
    best_x, best_f = copy(population[best_idx]), scores[best_idx]
    for generation in 1:max_generations
        order = sortperm(scores)
        next_population = [copy(population[order[1]])]
        while length(next_population) < population_size
            parent_a = population[_tournament(rng, scores)]
            parent_b = population[_tournament(rng, scores)]
            α = rand(rng, length(lb))
            child = α .* parent_a .+ (1 .- α) .* parent_b
            child .+= randn(rng, length(lb)) .* (ub .- lb) .* (0.08 * (1 - generation / max_generations))
            push!(next_population, _project_design(child, lb, ub))
        end
        population = next_population
        scores = [_objective(fitness_fn, x) for x in population]
        generation_best = argmin(scores)
        if scores[generation_best] < best_f
            best_x, best_f = copy(population[generation_best]), scores[generation_best]
        end
    end
    return OptResult(best_x, best_f, 1, (algorithm = :ga, generations = max_generations))
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
    swarm_size >= 2 && max_iterations >= 1 ||
        throw(ArgumentError("swarm_size must be ≥2 and iterations must be positive"))
    rng = MersenneTwister(seed)
    n = length(lb)
    positions = [_project_design(lb .+ rand(rng, n) .* (ub .- lb), lb, ub)
                 for _ in 1:swarm_size]
    velocities = [zeros(n) for _ in 1:swarm_size]
    scores = [_objective(fitness_fn, x) for x in positions]
    personal = deepcopy(positions)
    personal_scores = copy(scores)
    best_idx = argmin(scores)
    best_x, best_f = copy(positions[best_idx]), scores[best_idx]
    for iteration in 1:max_iterations
        inertia = 0.9 - 0.5 * iteration / max_iterations
        for i in 1:swarm_size
            velocities[i] .= inertia .* velocities[i] .+
                1.7 .* rand(rng, n) .* (personal[i] .- positions[i]) .+
                1.7 .* rand(rng, n) .* (best_x .- positions[i])
            positions[i] = _project_design(positions[i] .+ velocities[i], lb, ub)
            scores[i] = _objective(fitness_fn, positions[i])
            if scores[i] < personal_scores[i]
                personal[i], personal_scores[i] = copy(positions[i]), scores[i]
            end
            if scores[i] < best_f
                best_x, best_f = copy(positions[i]), scores[i]
            end
        end
    end
    return OptResult(best_x, best_f, 1, (algorithm = :pso, iterations = max_iterations))
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
    max_iterations >= 1 || throw(ArgumentError("max_iterations must be positive"))
    rng = MersenneTwister(seed)
    current = _project_design(copy(x0), lb, ub)
    current_f = _objective(fitness_fn, current)
    best, best_f = copy(current), current_f
    for iteration in 1:max_iterations
        temperature = max(1e-6, 1 - iteration / max_iterations)
        candidate = _project_design(current .+ randn(rng, length(lb)) .* (ub .- lb) .* temperature .* 0.1, lb, ub)
        candidate_f = _objective(fitness_fn, candidate)
        if candidate_f < current_f || rand(rng) < exp((current_f - candidate_f) / temperature)
            current, current_f = candidate, candidate_f
        end
        if current_f < best_f
            best, best_f = copy(current), current_f
        end
    end
    return OptResult(best, best_f, 1, (algorithm = :sa, iterations = max_iterations))
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
    max_iterations >= 1 || throw(ArgumentError("max_iterations must be positive"))
    x = _project_design(copy(x0), lb, ub)
    score = _objective(fitness_fn, x)
    step = (ub .- lb) ./ 4
    iterations = 0
    while iterations < max_iterations && maximum(step) > 1e-5
        improved = false
        for dimension in eachindex(x), direction in (-1.0, 1.0)
            candidate = copy(x)
            candidate[dimension] += direction * step[dimension]
            candidate = _project_design(candidate, lb, ub)
            candidate_score = _objective(fitness_fn, candidate)
            if candidate_score < score
                x, score, improved = candidate, candidate_score, true
            end
        end
        improved || (step ./= 2)
        iterations += 1
    end
    return OptResult(x, score, 1, (algorithm = :ps, iterations = iterations))
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
    length(lb) == length(ub) || throw(DimensionMismatch("lower and upper bounds must match"))
    all(lb .< ub) || throw(ArgumentError("each lower bound must be below its upper bound"))
    if optimizer == :ga
        return run_ga(fitness_fn; lb, ub, kwargs...)
    elseif optimizer == :pso
        return run_pso(fitness_fn; lb, ub, kwargs...)
    elseif optimizer == :sa
        x0 === nothing && throw(ArgumentError("x0 is required for simulated annealing"))
        return run_sa(fitness_fn, x0; lb, ub, kwargs...)
    elseif optimizer == :ps
        x0 === nothing && throw(ArgumentError("x0 is required for pattern search"))
        return run_ps(fitness_fn, x0; lb, ub, kwargs...)
    end
    throw(ArgumentError("unknown optimizer \"$optimizer\"; choose :ga, :pso, :sa, or :ps"))
end

function _project_design(x::AbstractVector, lb::Vector{Float64}, ub::Vector{Float64})
    length(x) == length(lb) == length(ub) ||
        throw(DimensionMismatch("design and bounds dimensions must match"))
    projected = clamp.(Float64.(x), lb, ub)
    if length(projected) >= 3
        projected[2] = min(projected[2], max(lb[2], 0.75 - projected[1]))
        projected[3] = min(projected[3], max(lb[3], 0.75 - projected[1]))
    end
    length(projected) >= 4 && (projected[4] = clamp(round(projected[4]), lb[4], ub[4]))
    return projected
end

function _objective(fitness_fn::Function, x::Vector{Float64})
    value = fitness_fn(x)
    value isa Number || throw(ArgumentError("optimizer objective must return a number"))
    return isfinite(value) ? Float64(value) : Inf
end

function _tournament(rng::AbstractRNG, scores::Vector{Float64})
    a, b = rand(rng, eachindex(scores)), rand(rng, eachindex(scores))
    return scores[a] <= scores[b] ? a : b
end
