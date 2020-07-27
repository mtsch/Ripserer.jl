module BenchTorus
using Ripserer

using BenchmarkTools
using Random
Random.seed!(1337)
include(joinpath(@__DIR__, "../test/data.jl"))
suite = BenchmarkGroup()

for (npoints, dim_max, threshold, modulus) in ((1024, 1, 1, 2),
                                               (1024, 1, nothing, 7),
                                               (256, 2, 1, 7),
                                               (256, 2, nothing, 2),
                                               (64, 4, nothing, 2),
                                               (81, 4, 1, 3))
    name = "n=$npoints, dim_max=$dim_max, modulus=$modulus, threshold=$threshold"
    if isnothing(threshold)
        bench = @benchmarkable ripserer($(rand_torus(npoints)),
                                        dim_max=$dim_max,
                                        modulus=$modulus)
    else
        bench = @benchmarkable ripserer($(rand_torus(npoints)),
                                        dim_max=$dim_max,
                                        modulus=$modulus,
                                        threshold=$threshold)
    end
    suite[name] = bench
end
end

BenchTorus.suite
