module BenchAlpha

using BenchmarkTools
using Ripserer

include(joinpath(@__DIR__, "utils.jl"))

benchmarks = [
    (file="annulus10k.pts", dim_max=2),
    (file="torus10k.pts", dim_max=2),
    (file="klein500.pts", dim_max=3),
]

suite = BenchmarkGroup()

for b in benchmarks
    name = b.file
    data = load_data(joinpath(@__DIR__, "data", b.file))
    suite[name] = @benchmarkable ripserer(
        Alpha($data); dim_max=$(b.dim_max)
    ) samples=1
end

end
BenchAlpha.suite
