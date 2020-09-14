module BenchCubical

using BenchmarkTools
using Ripserer

include(joinpath(@__DIR__, "utils.jl"))

benchmarks = [
    (file="bonsai64.dipha", dim_max=3),
    (file="bonsai128.dipha", dim_max=3),
    (file="lena1024.dipha", dim_max=2),
    (file="lena2048.dipha", dim_max=2),
]

suite = BenchmarkGroup()

for b in benchmarks
    name = b.file
    data = load_data(joinpath(@__DIR__, "data", b.file))
    suite[name] = @benchmarkable ripserer(
        Cubical($data); dim_max=$(b.dim_max)
    ) gcsample=true seconds=60
end

end
BenchCubical.suite
