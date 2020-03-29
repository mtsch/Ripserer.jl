using Ripserer
using BenchmarkTools
using SparseArrays

const SUITE = BenchmarkGroup()

dist = sparse(ones(1000, 1000))
for i in 1:size(dist, 1)
    dist[i, i] = 0
end
dropzeros!(dist)

col = Ripserer.CurrentColumn{3}()
binom = Ripserer.Binomial(1000, 1000)

SUITE["initialization"] = @benchmarkable Ripserer.initialize!(
    $col, $(DiameterSimplex{3}(1.0, 10, 1)), 5, $dist, $binom)
