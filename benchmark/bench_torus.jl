module BenchTorus
using Ripserer
using BenchmarkTools
include(joinpath(@__DIR__, "../test/data.jl"))
suite = BenchmarkGroup()

for (npoints, dim_max, threshold) in ((1024, 1, 1),
                                      (1024, 1, nothing),
                                      (256, 2, 1),
                                      (256, 2, nothing),
                                      (64, 4, nothing),
                                      (128, 4, 1))
    for modulus in (2, 7)
        name = "n=$npoints, dim_max=$dim_max, modulus=$modulus, threshold=$threshold"
        if isnothing(threshold)
            bench = @benchmarkable ripserer($(torus(npoints)),
                                            dim_max=$dim_max,
                                            modulus=$modulus)
        else
            bench = @benchmarkable ripserer($(torus(npoints)),
                                            dim_max=$dim_max,
                                            modulus=$modulus,
                                            threshold=$threshold)
        end
        suite[name] = bench
    end
end
end

BenchTorus.suite
