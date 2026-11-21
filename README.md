# Wing Structural Optimization

The supported package is **WingOpt**, a Julia implementation of the wing
structural-analysis and design-optimization workflow. The original MATLAB
scripts are retained as a historical reference; running the Julia package
does not require MATLAB.

## Requirements and installation

Install Julia 1.10 or later, then from the repository root run:

```sh
julia --project=WingOpt -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

The package and its project metadata live in `WingOpt/`. Airfoil coordinate
files used by the built-in case are in `foils/`.

## Quick start

```julia
using WingOpt

result, data = run_optimization(
    1;
    optimizer = :ps,
    max_iterations = 100,
)

result.x       # [front-web position, tip width, root width, wall mm, reinforcer m]
result.fval    # objective (minimized)
data.mass      # estimated complete-wing spar mass, kg
data.min_ms    # minimum calculated failure margin
```

Run Julia with the repository root as the working directory to use its `foils/`
directory, or pass `foil_dir` to `run_optimization`. The workflow also accepts
caller-provided `Vector{AeroLoad}` distributions. If no loads are supplied,
the entry point uses a documented example elliptic distribution with a 100 N
total lift; it is intended for smoke tests and examples, not as a substitute
for aerodynamic analysis.

Available optimizers are `:ga`, `:pso`, `:sa`, and `:ps`. The derivative-free
searches are implemented in the Julia package and use seeded random generators
for reproducibility. `:sa` and `:ps` accept the initial design through the
common optimizer interface; `run_optimization` supplies a feasible initial
design automatically.

## Package layout

`WingOpt/src/` contains typed configuration and result models and separate
modules for numerical utilities, materials, foil and section geometry, loads,
failure criteria, fitness, optimizers, and plots. `WingOpt/test/runtests.jl`
contains unit and workflow smoke tests. Run tests with `Pkg.test()` as shown
above.

The public entry point `run_optimization` builds run case 1, reads a Selig
`.dat` foil directly, constructs the reinforced-box spar, evaluates load
cases, and optimizes the design without `.mat` workspace files. Core routines
such as `parse_foil`, `generate_wing_stations`, `reinforced_box_section`,
`internal_loads_SS`, and the failure criteria are also available individually.

## Migration status and engineering limitations

The Julia implementation preserves the intended units and the documented
run-case-1 geometry, but MATLAB parity is not yet demonstrated against a
trusted reference result. In the archived sources, the fitness workflow
references `data/loads.mat` (not present in the repository), and the fixed
support (`FX`) implementation is incomplete. Therefore:

- The validated end-to-end path is run case 1 with a reinforced-box section
  and simply-supported (`ss`) loads.
- Fixed-support (`fx`) analysis fails explicitly rather than returning an
  invented result.
- Report-45 is not evaluated: the archived wing configuration lacks the
  design dive speed and aileron-span values that criterion requires.
- The default aerodynamic distribution is synthetic. Engineering results
  require user-supplied aerodynamic loads and independent validation.
- The MATLAB scripts are retained for reference and have not been deprecated.

Optimization results depend on the load case, material data, chosen optimizer,
and stopping settings. The example load is useful for exercising the package,
not for sizing an aircraft structure.

## References

- Greco, M. *Resistência dos Materiais: Uma Abordagem Sintética*.
- Megson, T. H. G. *Aircraft Structures for Engineering Students*.
- Airframe and Equipment Engineering Report No. 45, AD-A955 270, 1955.
