"""
    WingOpt

Julia translation of the Matlab wing structural optimization toolbox.

Original Matlab code: vsdsantos/wing-structural-optimization
Migration: Matlab → Julia (see migration plan in repository root).

## Quick start
```julia
using WingOpt
result, data = run_optimization(1; optimizer=:ga, max_generations=50)
```

## Module layout
| Source file            | Matlab equivalent                          |
|------------------------|--------------------------------------------|
| types.jl               | (all struct definitions)                   |
| utils.jl               | –                                          |
| material.jl            | material/material_properties.m             |
| geometry/transform.jl  | geometry/transform/*.m                     |
| geometry/section.jl    | geometry/section/*.m                       |
| geometry/foil.jl       | geometry/foil/{parse_foil,foil_properties}.m|
| geometry/stations.jl   | geometry/generate_wing_stations.m          |
| torsion.jl             | torsion/*.m + criterias/displacement/*.m   |
| sections.jl            | sections/*.m                               |
| parametrics.jl         | parametrics/*.m                            |
| forces.jl              | forces/internal_loads_SS.m                 |
| criteria/failure.jl    | criterias/failure/*.m                      |
| criteria/criteria.jl   | criterias/{safety_margin,beam_mass,report_45}.m|
| fitness.jl             | fit_function.m + weight_function.m         |
| config.jl              | run_cases.m + wing_geometry.m              |
| optimizers.jl          | optimizers/*.m                             |
| plots.jl               | plots/*.m                                  |
| main.jl                | start_opt.m                                |
"""
module WingOpt

using LinearAlgebra
using Interpolations
using Metaheuristics
using Optim

# ── Core types ────────────────────────────────────────────────────────────────
include("types.jl")

# ── Numerical utilities ───────────────────────────────────────────────────────
include("utils.jl")

# ── Physics sub-modules (order matters: each depends on the ones above) ───────
include("material.jl")

include("geometry/transform.jl")
include("geometry/section.jl")
include("geometry/foil.jl")
include("geometry/stations.jl")

include("torsion.jl")
include("sections.jl")
include("parametrics.jl")
include("forces.jl")

include("criteria/failure.jl")
include("criteria/criteria.jl")

include("fitness.jl")
include("config.jl")
include("optimizers.jl")

# ── Visualization (optional – loaded last so Plots is not required for headless use)
include("plots.jl")

# ── Entry point ───────────────────────────────────────────────────────────────
include("main.jl")

# ── Public API exports ────────────────────────────────────────────────────────
export
    # Types
    MaterialProperties,
    WingGeometry,
    RunConfig,
    FoilData,
    FoilGeometry,
    Station,
    SectionProperties,
    StructuredGeometry,
    AeroLoad,
    HalfLoads,
    InternalLoads,
    OptResult,
    FitData,

    # Material
    material_properties,

    # Geometry
    translate_inertia_tensor,
    rotate_inertia_tensor,
    principal_inertia_tensor,
    multielement_centroid,
    area_inertia_tensor,
    first_moment_of_area,
    section_centroid,
    parse_foil,
    foil_properties,
    generate_wing_stations,

    # Torsion / displacement
    integral_linha_media,
    ang_torcao_parede_fina,
    bicelular,
    wing_torsion,
    wing_displacement,

    # Sections
    reinforced_box_section,
    box_section,
    D_section,
    O_section,
    foil_shaped_section,
    compute_section,

    # Parametrics
    reinforced_box_parametric_geometry,

    # Forces
    internal_loads_SS,
    compute_internal_loads,

    # Failure criteria
    max_tension_criteria,
    tsai_wu_criteria,
    tsai_hill_criteria,

    # Structural criteria
    safety_margin,
    beam_mass,
    report_45,

    # Fitness
    weight_function,
    fit_function,
    example_aero_loads,

    # Config
    build_run_config,
    build_wing_geometry,

    # Optimisers
    run_ga,
    run_pso,
    run_sa,
    run_ps,
    run_optimizer,
    DEFAULT_LB,
    DEFAULT_UB,
    AINEQ,
    BINEQ,

    # Plots
    plot_internal_loads,
    plot_stiffness,
    plot_box_section,
    plot_spar_topview,

    # Entry point
    run_optimization

end # module WingOpt
