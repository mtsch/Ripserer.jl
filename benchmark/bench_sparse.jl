module BenchSparse
using Ripserer
using BenchmarkTools
using Random
include(joinpath(@__DIR__, "../test/data.jl"))
Random.seed!(666)
suite = BenchmarkGroup()

for modulus in (2, 17)
    suite["10 tori with 16 points, modulus=$modulus"] =
        @benchmarkable ripserer($(disconnected_tori(16, 10)),
                                dim_max=2, modulus=$modulus)
    suite["20 tori with 128 points, modulus=$modulus"] =
        @benchmarkable ripserer($(disconnected_tori(128, 20)),
                                dim_max=1, modulus=$modulus)
    suite["10 tori with 256 points, modulus=$modulus"] =
        @benchmarkable ripserer($(disconnected_tori(512, 10)),
                                dim_max=1, modulus=$modulus)
end

end

BenchSparse.suite
