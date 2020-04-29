module BenchSparse
using Ripserer
using BenchmarkTools
using Random
include(joinpath(@__DIR__, "../test/data.jl"))
Random.seed!(666)
suite = BenchmarkGroup()

for (n_points, n_tori, dim, modulus) in [(256, 20, 1, 2),
                                         (512, 10, 1, 3),
                                         (128, 10, 2, 3),
                                         (64, 20, 2, 2),
                                         (64, 10, 3, 2)]
    suite["$n_tori tori with $n_points points, dim_max=$dim, modulus=$modulus"] =
        @benchmarkable ripserer($(disconnected_tori(n_tori, n_points)),
                                dim_max=$dim, modulus=$modulus)
end
end

BenchSparse.suite
