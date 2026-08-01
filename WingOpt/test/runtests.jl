"""
Unit tests for WingOpt.jl
Run with:  julia --project=WingOpt WingOpt/test/runtests.jl
or:        cd WingOpt && julia --project -e 'using Pkg; Pkg.test()'
"""

using Test
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using WingOpt

@testset "WingOpt" begin
    @test_broken false  # placeholder: tests will be added as modules are implemented
end
