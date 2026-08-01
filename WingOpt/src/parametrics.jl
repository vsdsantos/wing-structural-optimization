"""
Parametric geometry generation for the spar cross-section.
Translates `parametrics/basic_box_section.m` and
`parametrics/reinforced_box_parametric_geometry.m`.
"""

"""
    reinforced_box_parametric_geometry(X, wing) -> StructuredGeometry

Build the spar parametric geometry from the optimisation variable vector `X`.

Variable mapping (1-based, same as Matlab):
    X[1] – web_pos_tip     (fraction of chord; same at all spans)
    X[2] – width_tip       (fraction of chord)
    X[3] – width_root      (fraction of chord; used at root and mid)
    X[4] – box_thickness   (integer 1–4; multiplied by a unit later in criteria)
    X[5] – reinforcer_thickness (m)

`wing.chord_transitions` = [root_chord, mid_chord, tip_chord]
`wing.span_transitions`  = [root_span,  mid_span,  tip_span ]

Returns a `StructuredGeometry` with one value per span break (3 control points).
"""
function reinforced_box_parametric_geometry(X::AbstractVector, wing::WingGeometry)
    error("Not implemented")
end
