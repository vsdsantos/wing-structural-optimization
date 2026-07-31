"""
Airfoil parsing and geometric property evaluation.
Replaces `geometry/foil/parse_foil.m` and `foil_properties.m`.

Key differences from the Matlab version:
- Uses Interpolations.jl LinearInterpolation instead of the Curve Fitting Toolbox.
- Returns in-memory FoilData instead of writing a .mat cache.
- `foil_properties` accepts a FoilData struct instead of a file name.
"""

"""
    parse_foil(foil_name; foil_dir="foils") -> FoilData

Parse a Selig-format `.dat` airfoil file and build linear interpolants for
the upper (extrados) and lower (intrados) surfaces.

The file is expected to live at `foil_dir/foil_name.dat`.
The first non-numeric line (title) is skipped automatically.
"""
function parse_foil(foil_name::String; foil_dir::String="foils")
    filename = joinpath(foil_dir, foil_name * ".dat")
    isfile(filename) || error("Airfoil file not found: $filename")

    # Read all numeric rows (skip header / comment lines)
    rows = Vector{NTuple{2,Float64}}()
    open(filename) do fh
        for line in eachline(fh)
            parts = split(strip(line))
            length(parts) < 2 && continue
            x = tryparse(Float64, parts[1])
            y = tryparse(Float64, parts[2])
            (x === nothing || y === nothing) && continue
            push!(rows, (x, y))
        end
    end
    isempty(rows) && error("No numeric data found in $filename")

    dat = Matrix{Float64}(undef, length(rows), 2)
    for (i, (x, y)) in enumerate(rows)
        dat[i, 1] = x;  dat[i, 2] = y
    end

    # Selig format: data runs TE→LE (upper) then LE→TE (lower).
    # Split at the leading edge (minimum x).
    _, i_le = findmin(dat[:, 1])

    ext_raw = dat[1:i_le, :]      # upper surface  TE→LE  (x decreasing)
    int_raw = dat[i_le:end, :]    # lower surface  LE→TE  (x increasing)

    # Reverse upper so x is monotonically increasing (required by Interpolations)
    ext_sorted = ext_raw[end:-1:1, :]
    int_sorted = int_raw

    # Build linear interpolants  x ∈ [0, 1] → y
    ex_lin = LinearInterpolation(ext_sorted[:, 1], ext_sorted[:, 2])
    in_lin = LinearInterpolation(int_sorted[:, 1], int_sorted[:, 2])

    # Evaluation step: 10 % of the minimum data spacing (matching Matlab)
    esp_ex = minimum(abs.(diff(ext_sorted[:, 1])))
    esp_in = minimum(abs.(diff(int_sorted[:, 1])))
    esp    = min(esp_ex, esp_in) * 1e-1

    return FoilData(ex_lin, in_lin, ext_sorted, int_sorted, esp)
end

# ---------------------------------------------------------------------------

"""
    foil_properties(chord, scale_thick, limits, offset, foil) -> FoilGeometry

Evaluate airfoil geometry at a spanwise station.

Arguments:
- `chord`       – local chord (m)
- `scale_thick` – thickness scaling offset (set to 0 for normal use)
- `limits`      – [a, b] front/back limits in fractional chord [0, 1]
- `offset`      – [upper_offset, lower_offset] y offsets (use [0,0] normally)
- `foil`        – FoilData from `parse_foil`

Returns a FoilGeometry with upper/lower surface arrays in the quarter-chord
coordinate system (x = 0 at c/4, positive toward TE; y = 0 at camber line).
"""
function foil_properties(chord::Float64, scale_thick::Float64,
                          limits::AbstractVector, offset::AbstractVector,
                          foil::FoilData)
    esp = foil.esp
    a, b = Float64(limits[1]), Float64(limits[2])

    # --- Camber-line shift: evaluate full chord to find y at c/4 ---
    xs_full = range(0.0, 1.0; step=esp)
    ext_full_y = foil.ex_lin.(xs_full)
    int_full_y = foil.in_lin.(xs_full)
    # scale & shift to quarter-chord frame
    x_full = collect(xs_full) .* chord .- chord / 4
    _, i_qc = findmin(abs.(x_full))        # index closest to x = 0
    shift = (ext_full_y[i_qc] + int_full_y[i_qc]) / 2 * chord

    # --- Spar section: evaluate at a→b ---
    xs = collect(range(a, b; step=esp))
    length(xs) < 2 && (xs = collect(LinRange(a, b, 10)))

    # Upper surface
    ext = hcat(xs, foil.ex_lin.(xs))
    ext = hcat(ext[:, 1] .* chord .- chord / 4,
               ext[:, 2] .* chord .- shift .- Float64(offset[1]))

    # Lower surface
    int = hcat(xs, foil.in_lin.(xs))
    int = hcat(int[:, 1] .* chord .- chord / 4,
               int[:, 2] .* chord .- shift .+ Float64(offset[2]))

    # Optional thickness scaling (kept for completeness; use scale_thick=0)
    if scale_thick != 0.0
        ext = _apply_scale_thick(ext,  scale_thick)
        int = _apply_scale_thick(int, -scale_thick)
    end

    # --- Arc lengths ---
    arc_len(pts) = sum(sqrt.(diff(pts[:, 1]).^2 .+ diff(pts[:, 2]).^2))
    ex_len    = arc_len(ext)
    in_len    = arc_len(int)
    front_len = ext[1, 2]   - int[1, 2]
    back_len  = ext[end, 2] - int[end, 2]

    area = trapz(int[:, 1], ext[:, 2] .- int[:, 2])

    x_all = vcat(ext[:, 1], int[:, 1])
    y_all = vcat(ext[:, 2], int[:, 2])

    return FoilGeometry(ext, int, ex_len, in_len, front_len, back_len,
                        ex_len + in_len + front_len + back_len, area, x_all, y_all)
end

# Internal: apply thickness scaling offset along the surface normal
function _apply_scale_thick(pts::Matrix{Float64}, scale::Float64)
    dx  = diff(pts[:, 1]);  dy  = diff(pts[:, 2])
    mid = pts[1:end-1, :] .+ hcat(dx, dy) ./ 2
    nrm = hcat(-dy, dx) ./ sqrt.(dx.^2 .+ dy.^2)   # outward normal
    return vcat(mid .+ nrm .* scale,
                reshape(pts[end, :], 1, 2))
end
