using WingOpt

semispan = 4.0
total_lift = 120.0
loads = WingOpt.Aerodynamics.elliptic_load(semispan, total_lift, 16)
panel_width = semispan / length(loads.x)
integrated_lift = 2 * sum(loads.dLdb) * panel_width

mktempdir() do dir
    result_path = joinpath(dir, "smoke_result.toml")
    result = Dict(
        "semispan_m" => semispan,
        "total_lift_n" => integrated_lift,
        "panel_count" => length(loads.x),
    )
    WingOpt.InputOutput.write_data(result_path, result)
    saved_result = WingOpt.InputOutput.read_data(result_path)
    @assert saved_result["total_lift_n"] ≈ total_lift
    println("WingOpt smoke test passed: ", saved_result)
end
