"""
Material property database.
Translates `material_properties.m` (switch/case → if/elseif).
"""

"""
    material_properties(name) -> MaterialProperties

Return mechanical properties for a named material.

Supported materials:
- `"carbono_bi"`  – bi-axial carbon fibre laminate
- `"carbono_uni"` – unidirectional carbon fibre laminate
- `"balsa"`       – balsa wood
- `"freijo"`      – freijo (hardwood) rod
"""
function material_properties(name::String)
    if name == "carbono_bi"
        return MaterialProperties(
            6.82e8,  3.41e8,  5.10e10, 0.28,   # axial
            6.82e8,  3.41e8,  4.8929e10, 0.02,  # transverse
            2.00e9,  7.20e7,  1545.0             # G, shear, density
        )
    elseif name == "carbono_uni"
        return MaterialProperties(
            1.13e9,  3.41e8,  9.43e10, 0.33,
            4.30e8,  1.33e8,  3.25e9,  0.03,
            1.98e9,  6.00e7,  1500.0
        )
    elseif name == "balsa"
        return MaterialProperties(
            4.00e7,  6.90e6,  4.10e9,  0.25,
            8.00e5,  1.00e5,  9.50e7,  0.25,
            1.66e8,  1.30e6,  150.0
        )
    elseif name == "freijo"
        return MaterialProperties(
            1.02e8,  4.50e7,  1.60e10, 0.30,
            4.20e6,  2.00e6,  9.00e8,  0.30,
            3.50e9,  8.70e6,  362.0
        )
    else
        error("Unknown material: \"$name\". " *
              "Valid names: carbono_bi, carbono_uni, balsa, freijo.")
    end
end
