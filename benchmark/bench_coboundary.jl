module BenchCoboundary
using Ripserer
using Ripserer: Coboundary
using BenchmarkTools
suite = BenchmarkGroup()

function count_cofaces(coboundary, sx)
    count = 0
    for i in 1:10000
        for coface in coboundary(sx, 2)
            count += 1
        end
    end
    count / 10000
end
dists = torus(4096)
sx = Simplex{2}(dists[1, 2], 1, 1)

coboundary_full = Coboundary(RipsFiltration(t), 2)
coboundary_sparse = Coboundary(SparseRipsFiltration(t), 2)

suite["full"] = @benchmarkable count_cofaces($coboundary_full, $sx)
suite["sparse"] = @benchmarkable count_cofaces($coboundary_sparse, $sx)

end

BenchTorus.suite
