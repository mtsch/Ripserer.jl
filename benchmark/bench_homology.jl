module BenchHomology

using BenchmarkTools
using Ripserer

include(joinpath(@__DIR__, "utils.jl"))

benchmarks = [
    (file="klein200.pts", dim_max=1, threshold=1.7, filtration=Rips),
    (file="sphere100.pts", dim_max=2, threshold=nothing, filtration=Rips),
    (file="lena1024.dipha", dim_max=2, threshold=nothing, filtration=Cubical),
    (file="bonsai64.dipha", dim_max=3, threshold=nothing, filtration=Cubical),
]

suite = BenchmarkGroup()

for b in benchmarks
    name = b.file
    data = load_data(joinpath(@__DIR__, "data", b.file))
    if isnothing(b.threshold)
        filtration = b.filtration(data)
    else
        filtration = b.filtration(data; threshold=b.threshold)
    end
    suite[name] = @benchmarkable(
        ripserer($filtration; dim_max=$(b.dim_max), alg=:involuted), samples = 1
    )
end

end
BenchHomology.suite
