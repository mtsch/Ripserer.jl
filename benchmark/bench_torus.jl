module BenchTorus
using Ripserer
using BenchmarkTools
include(joinpath(@__DIR__, "../test/data.jl"))
suite = BenchmarkGroup()

for (npoints, dim_max) in ((1024, 1), (128, 2), (64, 4))
    for modulus in (2, 7)
        for threshold in (nothing, 1)
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
end

BenchTorus.suite
