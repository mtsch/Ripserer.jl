module BenchProjectivePlane
using Ripserer
using BenchmarkTools
include(joinpath(@__DIR__, "../test/data.jl"))
suite = BenchmarkGroup()

suite["modulus=2"] =
    @benchmarkable ripserer(projective_plane, dim_max=2, modulus=2)
suite["modulus=3"] =
    @benchmarkable ripserer(projective_plane, dim_max=2, modulus=3)
suite["modulus=17"] =
    @benchmarkable ripserer(projective_plane, dim_max=2, modulus=17)

end

BenchProjectivePlane.suite
