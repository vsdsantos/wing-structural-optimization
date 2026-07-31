"""
Parametric geometry generation for the spar cross-section.
Translates `parametrics/basic_box_section.m` and
`parametrics/reinforced_box_parametric_geometry.m`.
"""

# ---------------------------------------------------------------------------
# Intermediate helper (parametric station)
# ---------------------------------------------------------------------------

"""
Internal helper that holds geometry data for one span break point
(used only during parametric construction).
"""
struct _ParamStation
    position::Float64   # front-web position (fraction of chord)
    width::Float64      # spar width (fraction of chord)
    thickness::Float64  # box wall thickness (m)
    reinforcer::Float64 # reinforcer thickness (m)
    chord::Float64      # local chord (m)
end

# ---------------------------------------------------------------------------
# basic_box_section
# ---------------------------------------------------------------------------

"""
    basic_box_section(root, tip) -> NamedTuple

Interpolate box-section geometry between two span-break stations.

Returns a NamedTuple with fields:
- `x_position` – [root, tip] front-web positions (m, offset from c/4)
- `width`       – [root, tip] spar widths (m)
- `thickness`   – [root, tip] wall thicknesses (m)
- `reinforcer`  – [root, tip] reinforcer thicknesses (m)

(Replicates the Matlab `section.limits` field but drops the cell-array
wrapping; `generate_wing_stations` reconstructs limits from x_position + width.)
"""
function basic_box_section(root::_ParamStation, tip::_ParamStation)
    x_position = [root.position * root.chord - root.chord / 4,
                   tip.position  * tip.chord  - tip.chord  / 4]
    width      = [root.width * root.chord,
                   tip.width  * tip.chord ]
    thickness  = [root.thickness, tip.thickness]
    reinforcer = [root.reinforcer, tip.reinforcer]
    return (x_position = x_position, width = width,
            thickness  = thickness,  reinforcer = reinforcer)
end

# ---------------------------------------------------------------------------
# reinforced_box_parametric_geometry
# ---------------------------------------------------------------------------

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
    # Round the integer variable
    box_thick  = round(Int, clamp(X[4], 1, 4)) * 1e-3   # integer ×1 mm in metres
    refor_thick = X[5]

    c_root = wing.chord_transitions[1]
    c_mid  = wing.chord_transitions[2]
    c_tip  = wing.chord_transitions[3]

    root = _ParamStation(X[1], X[3], box_thick, refor_thick, c_root)
    mid  = _ParamStation(X[1], X[3], box_thick, refor_thick, c_mid)
    tip  = _ParamStation(X[1], X[2], box_thick, refor_thick, c_tip)

    sec1 = basic_box_section(root, mid)
    sec2 = basic_box_section(mid,  tip)

    return StructuredGeometry(
        [sec1.x_position..., sec2.x_position[2]],
        [sec1.width...,      sec2.width[2]      ],
        [sec1.thickness...,  sec2.thickness[2]  ],
        [root.reinforcer, mid.reinforcer, tip.reinforcer]
    )
end
