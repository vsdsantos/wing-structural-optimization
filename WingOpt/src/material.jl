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
    error("Not implemented")
end
