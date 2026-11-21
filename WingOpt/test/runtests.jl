"""
Unit tests for WingOpt.jl
Run with:  julia --project=WingOpt WingOpt/test/runtests.jl
or:        cd WingOpt && julia --project -e 'using Pkg; Pkg.test()'
"""

using Test
using WingOpt

@testset "WingOpt scaffold" begin
    @test WingOpt.Geometry isa Module
    @test WingOpt.Loads isa Module
    @test WingOpt.StructuralAnalysis isa Module
    @test WingOpt.Aerodynamics isa Module
    @test WingOpt.Optimization isa Module
    @test WingOpt.InputOutput isa Module
    @test WingOpt.Visualization isa Module

    @testset "elliptic load" begin
        semispan = 4.0
        total_lift = 120.0
        panels = 16
        loads = WingOpt.Aerodynamics.elliptic_load(semispan, total_lift, panels)

        @test length(loads.x) == panels
        @test all(diff(loads.x) .> 0)
        @test 2 * sum(loads.dLdb) * semispan / panels ≈ total_lift
        @test all(iszero, loads.dMadb)
        @test_throws ArgumentError WingOpt.Aerodynamics.elliptic_load(0, 1, panels)
    end

    @testset "TOML data round trip" begin
        mktempdir() do dir
            path = joinpath(dir, "result.toml")
            data = Dict("lift_n" => 120.0, "panel_count" => 16)

            @test WingOpt.InputOutput.write_data(path, data) == path
            @test WingOpt.InputOutput.read_data(path) == data
        end
    end
end
