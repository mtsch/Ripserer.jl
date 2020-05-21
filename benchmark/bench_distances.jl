module BenchDistances
using Ripserer
using BenchmarkTools
using FileIO
suite = BenchmarkGroup()

image = Float32.(getfield.(load(joinpath(@__DIR__, "helix_nebula1200px.jpg")), :r))
_, pd_big = ripserer(Cubical(image[1:200, 1:200]))
_, pd_med = ripserer(Cubical(image[301:400, 1:200]))
_, pd_sml = ripserer(Cubical(image[301:400, 150:200]))
_, pd_min = ripserer(Cubical(image[351:400, 150:200]))

large = (pd_big, pd_med)
medium = (pd_med, pd_sml)
small = (pd_med, pd_min)

suite["Bottleneck() large"] = @benchmarkable matching(Bottleneck(), large[1], large[2])
suite["Bottleneck() medium"] = @benchmarkable matching(Bottleneck(), medium[1], medium[2])
suite["Bottleneck() small"] = @benchmarkable matching(Bottleneck(), small[1], small[2])
suite["Wasserstein() medium"] = @benchmarkable matching(Wasserstein(), medium[1], medium[2])
suite["Wasserstein() small"] = @benchmarkable matching(Wasserstein(), small[1], small[2])

end

BenchDistances.suite
