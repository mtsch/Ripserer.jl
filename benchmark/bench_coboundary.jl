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
# Distances are between 0 and 2.
dists = rand_dist_matrix(4000)
sx = Simplex{2}(dists[1, 2], 1, 1)

coboundary_full_nothreshold = Coboundary(RipsFiltration(dists, threshold=10), 2)
coboundary_full_threshold1 = Coboundary(RipsFiltration(dists, threshold=1), 2)
coboundary_sparse_75 = Coboundary(SparseRipsFiltration(dists, threshold=1.5), 2)
coboundary_sparse_50 = Coboundary(SparseRipsFiltration(dists, threshold=1), 2)
coboundary_sparse_25 = Coboundary(SparseRipsFiltration(dists, threshold=0.5), 2)

suite["full, no threshold"] =
    @benchmarkable count_cofaces($coboundary_full_nothreshold, $sx)
suite["full, threshold=1"] =
    @benchmarkable count_cofaces($coboundary_full_threshold1, $sx)
suite["sparse, 75% full"] =
    @benchmarkable count_cofaces($coboundary_sparse_75, $sx)
suite["sparse, 50% full"] =
    @benchmarkable count_cofaces($coboundary_sparse_50, $sx)
suite["sparse, 25% full"] =
    @benchmarkable count_cofaces($coboundary_sparse_25, $sx)

end

BenchCoboundary.suite
