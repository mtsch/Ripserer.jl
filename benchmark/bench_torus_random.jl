module BenchTorusRandom
using Ripserer
using BenchmarkTools
using Random
include(joinpath(@__DIR__, "../test/data.jl"))
suite = BenchmarkGroup()

Random.seed!(1337)

for (npoints, dim_max) in ((1024, 1), (128, 2), (64, 4))
    name = "n=$npoints, dim_max=$dim_max"
    suite[name] =
        @benchmarkable ripserer($(rand_torus(npoints)), dim_max=$dim_max)
end
end

BenchTorusRandom.suite
