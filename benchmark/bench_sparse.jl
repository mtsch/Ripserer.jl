module BenchSparse
using Ripserer
using BenchmarkTools
using Random
include(joinpath(@__DIR__, "../test/filtrations/test-datasets.jl"))
Random.seed!(666)
suite = BenchmarkGroup()

for (n_points, n_tori, dim, modulus) in [(200, 30, 1, 2),
                                         (300, 20, 1, 3),
                                         (100, 20, 2, 3),
                                         (50, 20, 3, 2)]
    suite["$n_tori tori with $n_points points, dim_max=$dim, modulus=$modulus"] =
        @benchmarkable ripserer($(disconnected_tori(n_points, n_tori)),
                                dim_max=$dim, modulus=$modulus)
end
end

BenchSparse.suite
