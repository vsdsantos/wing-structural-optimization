# MATLAB workflow inventory and migration specification

This document records the MATLAB source and the intended Julia migration
contract. It distinguishes behavior stated in the MATLAB files from decisions
that still need confirmation. The MATLAB program is a collection of scripts
and functions rather than a self-contained, currently reproducible application.

## Entry points and intended workflow

| MATLAB entry point | Interface and responsibility | Outputs / next step |
|---|---|---|
| `start_opt.m` | Workspace script. Selects run case 1, prepares configuration and airfoil/wing data, sets optimizer to `ga` and `parallel=false`, then starts optimization. | Calls `run_cases.m`, `geometry/foil/parse_foil.m`, `wing_geometry.m`, an optimizer, and finally `fit_function(optresult.x)`. |
| `run_cases.m` | Workspace script; `run` selects a configured case. Only case 1 is defined. | Saves `data/run_config.mat`. |
| `geometry/foil/parse_foil.m` | Workspace script; reads the selected airfoil `.dat` and fits interpolation curves. | Saves `foils/<foil_name>.mat`. The comment in `start_opt.m` claiming it generates `data/foil_fit.mat` is stale. |
| `wing_geometry.m` | Workspace script; reads run configuration and builds the tapered, symmetric planform. | Saves `data/wing_geometry.mat`. |
| `fit_function.m` | Intended objective/evaluation function, `fit_function(X) -> [raw_fit, data]`. | Builds geometry, computes structural response, mass and fitness; returns detailed analysis data. |

Intended execution order:

1. `start_opt.m` chooses a run case and calls the configuration, foil parsing,
   and wing-geometry scripts.
2. The selected optimizer calls `fit_function(X)` for each candidate design.
3. `fit_function` is intended to create parametric geometry and span stations,
   calculate section properties, compute internal loads, evaluate safety
   margins and the aeroelastic report, integrate spar mass, and score the
   result using `weight_function.m`.
4. `start_opt.m` reevaluates the best design to populate `data`.

The scripts depend on the MATLAB current working directory, `addpath`, shared
workspace variables, and relative paths. The equivalent Julia interface is
intended to be an ordinary function without workspace or filesystem side
effects:

```julia
run_optimization(run_case=1;
    optimizer=:ga, parallel=false, foil_dir="foils",
    loads=nothing, n_panels=20, kwargs...) -> (OptResult, FitData)
```

This signature is documented in `WingOpt/src/main.jl`; its implementation is
currently a stub. Supported Julia optimizer names are `:ga`, `:pso`, `:sa`,
and `:ps`; `:pso` corresponds to MATLAB's `psw_opt.m`.

## MATLAB source inventory and Julia ownership

Every MATLAB source is assigned to its planned Julia owner. A dash means the
file is a script or has no independent planned module; its behavior belongs
to the listed owner.

| MATLAB source(s) | Responsibility | Planned Julia owner |
|---|---|---|
| `start_opt.m`, `run_cases.m`, `wing_geometry.m` | Workflow setup, run cases, wing configuration | `WingOpt/src/main.jl`, `config.jl` |
| `fit_function.m`, `weight_function.m` | Full evaluation pipeline and scalar objective | `fitness.jl` |
| `material/material_properties.m` | Material properties | `material.jl` |
| `geometry/foil/parse_foil.m`, `geometry/foil/foil_properties.m` | Airfoil parsing and foil-derived geometry | `geometry/foil.jl` |
| `geometry/generate_wing_stations.m` | Spanwise station interpolation | `geometry/stations.jl` |
| `geometry/section/area_inertia_tensor.m`, `geometry/section/first_moment_of_area.m`, `geometry/section/multielement_centroid.m`, `geometry/section/section_centroid.m` | Section integrals and centroids | `geometry/section.jl` |
| `geometry/transform/principal_inertia_tensor.m`, `geometry/transform/rotate_inertia_tensor.m`, `geometry/transform/translate_inertia_tensor.m` | Inertia tensor transformations | `geometry/transform.jl` |
| `parametrics/basic_box_section.m`, `parametrics/reinforced_box_parametric_geometry.m` | Design-variable-to-spar geometry | `parametrics.jl` |
| `sections/D_section.m`, `sections/O_section.m`, `sections/box_section.m`, `sections/foil_shaped_section.m`, `sections/reinforced_box_section.m` | Cross-section property calculations | `sections.jl` |
| `forces/internal_loads_FX.m`, `forces/internal_loads_SS.m`, `forces/panel_cumtrapz.m`, `forces/static_aeroelastic_loads.m` | Internal loads, integration, and aeroelastic load adjustment | `forces.jl` (`panel_cumtrapz` may share numerical utilities) |
| `torsion/ang_torcao_macico.m`, `torsion/ang_torcao_parede_fina.m`, `torsion/bicelular.m`, `torsion/integral_linha_media.m` | Torsional response and integrals | `torsion.jl` |
| `criterias/displacement/wing_displacement.m`, `criterias/displacement/wing_torsion.m` | Wing displacement and twist | `torsion.jl` |
| `criterias/beam_mass.m`, `criterias/safety_margin.m`, `criterias/aeroelastic/report_45.m`, `criterias/failure/max_tension_criteria.m`, `criterias/failure/tsai_hill_criteria.m`, `criterias/failure/tsai_wu_criteria.m` | Mass, failure/safety criteria, aeroelastic report | `criteria/criteria.jl`, `criteria/failure.jl` |
| `optimizers/ga_opt.m`, `optimizers/ps_opt.m`, `optimizers/psw_opt.m`, `optimizers/sa_opt.m` | Genetic algorithm, pattern search, particle swarm, simulated annealing | `optimizers.jl` |
| `plots/model_fem.m`, `plots/plot_esforcos.m`, `plots/plot_longarina.m`, `plots/plot_longarina_3d.m`, `plots/plot_prop_geo.m`, `plots/plot_sec_perfilada.m`, `plots/plot_sec_retangular.m`, `plots/stlwrite.m` | Analysis visualization, geometry rendering and STL export | `plots.jl` |

The table covers all 51 tracked `.m` files. Source-to-module mapping is a
migration target, not a claim that the Julia functions are implemented.

## Data assets

| Path | Role / status |
|---|---|
| `data/run_config.mat` | Saved output of the run-case setup. |
| `data/foil_fit.mat` | Tracked MATLAB artifact with no producer referenced by the current source; do not assume it is the active foil-fit fixture. |
| `data/wing_geometry.mat` | Saved output of `wing_geometry.m`. |
| `run/iter.mat` | Saved iteration artifact; its producing optimizer/run metadata is not recorded. |
| `foils/circle.dat`, `foils/s1210.dat`, `foils/s1223.dat`, `foils/tempest.dat` | Airfoil coordinate input files; case 1 selects `s1223`. |
| `data/loads.mat` | **Missing.** `fit_function.m` attempts to load this file and expects a `loads` array. |

The tracked MATLAB-saved artifacts can be used as partial migration references
where their producer is identifiable (`run_config.mat` and
`wing_geometry.mat`). They are not a complete analysis fixture: the repository
does not record their MATLAB release, source revision, or a set of expected
numeric values, and no aerodynamic load fixture is present.

## Numerical and API contract

### Case 1 and planform

`run_cases.m` defines the only active case:

| Setting | Case 1 |
|---|---|
| `section_type` | `reinforced_box` |
| `constraint_type` | `ss` (simply supported) |
| `constraint_gap` | `0.1` |
| `constraint_angle` | `0` |
| `foil_name` | `s1223` |

`wing_geometry.m` sets root chord `1`, full span `7`, taper ratios `[1, 0.5]`,
and taper position `0.75` of the semispan. The resulting span transitions are
`[0, 2.625, 3.5]`; chord transitions are `[1, 1, 0.5]`. Sweep, dihedral, and
geometric torsion are set to zero and marked not implemented.

The intended physical convention is a symmetric full-span wing represented
using semispan stations from the symmetry plane (`0`) to the tip (`3.5`).
The force code reflects stations for left/right-wing calculations. In its
comment, `xl` is measured relative to the quarter chord, positive toward the
trailing edge. `wing_geometry.m` says positive geometric torsion is pitch-up.
The exact axes/origin and handedness should be made explicit in the Julia
types before porting signed-load and tensor formulas.

The geometry constants imply metres, but the MATLAB sources do not declare a
global units schema. The migration should use SI internally and document
conversion at every external interface. Aerodynamic load units are not
specified in the repository.

### Optimization vector, bounds, and constraints

`start_opt.m` and `parametrics/reinforced_box_parametric_geometry.m` define
five design variables:

| Index | Meaning | Lower | Upper | Source convention |
|---:|---|---:|---:|---|
| 1 | Web/front position as chord fraction | 0.01 | 0.60 | Dimensionless |
| 2 | Tip box width as chord fraction | 0.01 | 0.40 | Dimensionless |
| 3 | Root/mid box width as chord fraction | 0.01 | 0.40 | Dimensionless |
| 4 | Box thickness | 1 | 4 | Integer variable (`ga` integer index `[4]`); comments imply 1 mm increments, but MATLAB code shows no unit conversion. |
| 5 | Reinforcer thickness | 0.002 | 0.010 | Julia draft labels metres; MATLAB source has no explicit unit declaration. |

Linear inequalities:

```text
X[1] + X[2] <= 0.75
X[1] + X[3] <= 0.75
```

MATLAB `ga_opt.m` uses 100 generations, population size 200, and linear
feasible creation. The optimizer is selected as `ga` and parallel evaluation
is disabled in `start_opt.m`. `ps_opt.m` and `sa_opt.m` set `x0=[]`, so their
current entry-point path does not provide a usable initial design. Other
solver defaults are MATLAB/toolbox defaults and are not pinned in source.

The Julia migration must resolve and document the box/reinforcer thickness
units before treating all design variables as SI. Preserve the integer
restriction on the fourth design variable; do not silently round continuous
optimizer output without a specified rule.

### Loads, output, objective, and errors

`fit_function.m` expects each `loads` entry to provide at least `x`, `dLdb`,
and `dMadb`; code also references `n_sec`. `internal_loads_SS.m` applies a
factor of `1.5 * 1.15` to distributed lift and torsion. This is part of the
legacy calculation and must remain explicit in the migrated contract rather
than being hidden in a fixture.

The MATLAB optimizer returns an `optresult` structure containing at least
`x`, `fval`, `exitflag`, and `output`; GA also returns `population` and
`score`. The intended Julia result is `OptResult` plus `FitData`, with
intermediate data including mass, margins, loads, section properties, and
the governing load case. Confirm exact Julia field names against
`WingOpt/src/types.jl` during implementation.

`weight_function.m` maximizes a gain by returning its negative as the
minimization objective. Its branches are:

| Condition | Gain before negation |
|---|---:|
| Safety margin `<= 0` | `-mass / (1 + safety_margin)` |
| `0 < safety_margin <= 1` and report `< 1` | `1.5` |
| Report `>= 1` | `0.5` |
| Otherwise | `1` |

The final score is `-gain / mass`. The first branch and boundary cases can
produce singular or counterintuitive values; preserve only after comparison
with a valid MATLAB baseline and explicitly define behavior for zero mass,
margin `-1`, non-finite loads, and other invalid inputs.

The current MATLAB files do not provide a reliable uniform failure contract.
The migration API should reject malformed dimensions, unknown run cases or
section/optimizer names, missing load fields, and non-finite inputs with
descriptive errors rather than relying on workspace-variable or indexing
failures.

## MATLAB and Julia dependencies

| Dependency | Used by | Migration note |
|---|---|---|
| MATLAB base language | Scripts, arrays, `load`/`save`, interpolation and numerical integration (`interp1`, `trapz`, `cumtrapz`, `linsolve`) | Reimplement deterministic math and explicit file/config handling in Julia. |
| Curve Fitting Toolbox | `fittype`, `fit`, and `prepareCurveData` in `geometry/foil/parse_foil.m` | Used for linear interpolation fits of airfoil coordinates; capture the interpolation behavior and spacing in a fixture before replacing it. |
| Global Optimization Toolbox | `ga`, `particleswarm`, `patternsearch`, `simulannealbnd`, `optimoptions`, solver callbacks and plotting functions in `optimizers/*.m` | Required by all four legacy optimizer wrappers. Record MATLAB/toolbox release, random seed, and solver settings for comparable optimizer baselines. |
| Parallel Computing Toolbox | Potentially required if the `UseParallel` options are enabled | `start_opt.m` sets `parallel=false`; parallel mode is not the default workflow. |
| Repository plotting/export code | `plots/*.m`, including `stlwrite.m` | Visualization/export is secondary to numeric parity and should not affect headless analysis. |
| Julia package dependencies | `WingOpt/Project.toml`: `Interpolations`, `Metaheuristics`, `Optim`, `Plots`, and standard-library `LinearAlgebra` | The corresponding planned owners are the Julia geometry, optimizer, and plotting modules. |

## Coupling, duplication, and known gaps

- The top-level process uses mutable workspace state, relative `load`/`save`
  paths, and `addpath`; config, loads, geometry, and optimizer output are
  coupled through files rather than function arguments.
- `fit_function.m` combines geometry, sections, loads, failure criteria,
  reporting, mass, objective scoring, and drawing data. Split these into the
  planned Julia owners and pass typed inputs/results between them.
- `fit_function.m` references `data/loads.mat`, which is absent, and symbols
  not defined in the tracked sources, including `n_sec`, `beam_forces_SS`,
  `beam_forces_FX`, `long_mass`, `station`, `sections`, `long_geo`, `freijo`,
  `balsa`, and `desenho`. It also calls `foild_shaped_section` (while the
  tracked function is `foil_shaped_section`).
- `fit_function.m`'s `beam_forces_*` calls do not match the tracked
  `internal_loads_*` names, and it switches over section types that its
  parametric geometry setup does not implement (only `reinforced_box` is
  handled there).
- `forces/internal_loads_FX.m` contains a malformed expression at
  `cumtrapz(s_0(,w)`, so the fixed-support path cannot run as checked in.
- `run_cases.m` configures only case 1; it has no default for unknown cases.
  `wing_geometry.m` explicitly marks sweep, dihedral, geometric torsion, and
  foil variation as unimplemented.
- Optimizer wrappers persist iteration state and configure many plotting
  callbacks; callback behavior, randomness, and solver defaults are tightly
  coupled to MATLAB toolbox releases.

These are migration blockers, not behaviors to silently mask. Fix or declare
out of scope each gap before claiming parity for that workflow.

## Baseline cases, fixtures, and tolerances

| Case | Workflow covered | Reference material and status |
|---|---|---|
| `CFG-01` | Run case and planform generation | Case 1 as specified above. `data/run_config.mat` and `data/wing_geometry.mat` are tracked saved MATLAB outputs and are partial candidate fixtures. Record MATLAB release and extract named field values before using them as golden numeric assertions. |
| `FOIL-01` | Airfoil parsing/fit | Parse `foils/s1223.dat` with `geometry/foil/parse_foil.m`; the expected output is `foils/s1223.mat`, which is not tracked. `data/foil_fit.mat` is not emitted by the current parser and is not accepted as its reference fixture without provenance. Capture fitted arrays and parsing settings from MATLAB before asserting numeric parity. |
| `ANALYSIS-01` | Section, loads, criteria, mass, and fitness | Intended case 1 using `X=[0.1, 0.2, 0.2, 2, 0.004]` is a proposed representative design, not an observed result. **No end-to-end reference output is available:** `data/loads.mat` is absent and `fit_function.m` has undefined dependencies. Keep this numeric case blocked until valid loads and missing implementation inputs are supplied. |
| `OPT-01` | Optimizer wrappers | GA case 1 uses the settings above. **No reproducible optimizer baseline is available:** no input loads, MATLAB/toolbox release, random seed, or complete objective execution is recorded. Optimizer parity is out of scope until `ANALYSIS-01` can run. |
| `PLOT-01` | Plotting and STL export | No numeric fixture is required; visual/export behavior is secondary and has no documented expected output. Treat as out of scope for numeric parity, but smoke-test it when implemented. |

The tracked `.mat` files are retained as partial reference assets, not newly
generated fixtures. A MATLAB or compatible runtime is not installed in the
environment used to prepare this inventory, so no MATLAB outputs were
regenerated or independently verified. `data/loads.mat` is also absent.
Consequently, a complete end-to-end MATLAB reference fixture remains an
explicit open deliverable; do not substitute fabricated values or Julia stub
outputs for it.

Proposed numeric comparison policy once named outputs are captured:

- Compare exact for categorical fields and integer design variables.
- For deterministic finite SI numeric outputs, begin with `rtol=1e-8` and
  `atol=1e-10`; document per-field units and justify any looser tolerance.
- For optimization, record MATLAB release, toolbox release, random seed,
  options, and initial design. Compare feasibility and objective with explicit
  tolerances; do not require identical optimizer trajectories or designs
  across solver versions.
- Each fixture should record source revision, MATLAB/toolbox versions, run
  case, optimizer settings, load input identity, and expected output fields.
