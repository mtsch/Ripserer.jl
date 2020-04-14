module BenchSparse
using Ripserer
using BenchmarkTools
using Random
include(joinpath(@__DIR__, "../test/data.jl"))
Random.seed!(666)
suite = BenchmarkGroup()

for modulus in (2, 17)
    for (n_points, n_tori, dim) in [(128,20,1), (256,10,1), (64,10,2), (32,20,2), (32,5,3)]
        suite["$n_tori tori with $n_points points, dim_max=$dim, modulus=$modulus"] =
            @benchmarkable ripserer($(disconnected_tori(n_tori, n_points)),
                                    dim_max=$dim, modulus=$modulus)
    end
end

end

BenchSparse.suite
