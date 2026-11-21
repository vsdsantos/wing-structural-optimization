module Geometry
using ..WingOpt: translate_inertia_tensor, rotate_inertia_tensor,
    principal_inertia_tensor, multielement_centroid, area_inertia_tensor,
    first_moment_of_area, section_centroid, parse_foil, foil_properties,
    generate_wing_stations
export translate_inertia_tensor, rotate_inertia_tensor, principal_inertia_tensor,
    multielement_centroid, area_inertia_tensor, first_moment_of_area,
    section_centroid, parse_foil, foil_properties, generate_wing_stations
end

module Loads
using ..WingOpt: AeroLoad, HalfLoads, InternalLoads, internal_loads_SS,
    compute_internal_loads
export AeroLoad, HalfLoads, InternalLoads, internal_loads_SS,
    compute_internal_loads
end

module StructuralAnalysis
using ..WingOpt: MaterialProperties, SectionProperties, safety_margin,
    beam_mass, report_45, max_tension_criteria, tsai_wu_criteria,
    tsai_hill_criteria, integral_linha_media, ang_torcao_parede_fina,
    bicelular, wing_torsion, wing_displacement
export MaterialProperties, SectionProperties, safety_margin, beam_mass,
    report_45, max_tension_criteria, tsai_wu_criteria, tsai_hill_criteria,
    integral_linha_media, ang_torcao_parede_fina, bicelular, wing_torsion,
    wing_displacement
end

module Optimization
using ..WingOpt: FitData, OptResult, fit_function, weight_function,
    run_ga, run_pso, run_sa, run_ps, run_optimizer, run_optimization
export FitData, OptResult, fit_function, weight_function, run_ga, run_pso,
    run_sa, run_ps, run_optimizer, run_optimization
end

module Visualization
using ..WingOpt: plot_internal_loads, plot_stiffness, plot_box_section,
    plot_spar_topview
export plot_internal_loads, plot_stiffness, plot_box_section, plot_spar_topview
end
