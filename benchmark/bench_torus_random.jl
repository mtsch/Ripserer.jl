module BenchTorusRandom
using Ripserer
using BenchmarkTools
using Random
include(joinpath(@__DIR__, "../test/data.jl"))
suite = BenchmarkGroup()

Random.seed!(1337)

for (npoints, dim_max) in ((1024, 1), (128, 2), (64, 4))
    for modulus in (2, 7)
        for threshold in (nothing, 1)
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
end
end

BenchTorusRandom.suite
