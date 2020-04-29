module BenchCoboundary
using Random

using Ripserer
using Ripserer: coboundary
using BenchmarkTools
suite = BenchmarkGroup()
include(joinpath(@__DIR__, "../test/data.jl"))
Random.seed!(7350)

function count_cofaces(filtration, sx; reps=10000)
    count = 0
    for i in 1:reps
        for coface in coboundary(filtration, sx)
            count += 1
        end
    end
    count ÷ reps
end
# Distances are between 0 and 2.
dists = rand_dist_matrix(4000)
dists[1, 2] = dists[2, 1] = 0.5
sx = Simplex{2, 2}(dists[1, 2], 1, 1)

flt_full_nothreshold = RipsFiltration(dists, threshold=10)
flt_full_threshold1 = RipsFiltration(dists, threshold=1)
flt_sparse_75 = SparseRipsFiltration(dists, threshold=1.5)
flt_sparse_50 = SparseRipsFiltration(dists, threshold=1)
flt_sparse_25 = SparseRipsFiltration(dists, threshold=0.5)

suite["full, no threshold"] =
    @benchmarkable count_cofaces($flt_full_nothreshold, $sx, reps=50_000)
suite["full, threshold=1"] =
    @benchmarkable count_cofaces($flt_full_threshold1, $sx, reps=50_000)
suite["sparse, 75% full"] =
    @benchmarkable count_cofaces($flt_sparse_75, $sx)
suite["sparse, 50% full"] =
    @benchmarkable count_cofaces($flt_sparse_50, $sx, reps=20_000)
suite["sparse, 25% full"] =
    @benchmarkable count_cofaces($flt_sparse_25, $sx, reps=40_000)

end

BenchCoboundary.suite
