module BenchRips

using BenchmarkTools
using Ripserer

include(joinpath(@__DIR__, "utils.jl"))

full_benchmarks = [
    (file="o3_1024.pts", dim_max=3, threshold=1.8),
    (file="fract-r.dist", dim_max=2, threshold=nothing),
    (file="dragon2000.pts", dim_max=1, threshold=nothing),
]
sparse_benchmarks = [(file="alpha_4_sphere_2000.spdist", dim_max=4)]

suite = BenchmarkGroup()
suite["sparse"] = BenchmarkGroup()
suite["dense"] = BenchmarkGroup()

for b in full_benchmarks
    for sparse in (true, false)
        name = b.file
        data = load_data(joinpath(@__DIR__, "data", b.file))
        suite[sparse ? "sparse" : "dense"][name] = @benchmarkable ripserer(
            $data; dim_max=$(b.dim_max), threshold=$(b.threshold), sparse=$sparse
        ) samples = 1
    end
end

for b in sparse_benchmarks
    name = b.file
    data = load_data(joinpath(@__DIR__, "data", b.file))
    suite["sparse"][name] = @benchmarkable ripserer($data; dim_max=$(b.dim_max)) samples = 2
end

end
BenchRips.suite
