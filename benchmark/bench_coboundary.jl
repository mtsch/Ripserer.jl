module BenchCoboundary
using Random

using Ripserer
using Ripserer: Coboundary
using BenchmarkTools
suite = BenchmarkGroup()
include(joinpath(@__DIR__, "../test/data.jl"))
Random.seed!(7350)

function count_cofaces(coboundary, sx)
    count = 0
    for i in 1:10000
        for coface in coboundary(sx, 2)
            count += 1
        end
    end
    count / 10000
end
dists = rand_dist_matrix(4000)
sx = Simplex{2}(dists[1, 2], 1, 1)

coboundary_full = Coboundary(RipsFiltration(dists, threshold=1), 2)
coboundary_sparse = Coboundary(SparseRipsFiltration(dists, threshold=1), 2)

suite["full"] =
    @benchmarkable count_cofaces($coboundary_full, $sx)
suite["sparse"] =
    @benchmarkable count_cofaces($coboundary_sparse, $sx)

end

BenchTorus.suite
