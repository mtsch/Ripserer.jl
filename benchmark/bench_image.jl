module BenchImage
using BenchmarkTools
using FileIO
using Ripserer
hubble = getfield.(load(joinpath(@__DIR__, "hubble1024px.jpg")), :r)
helix = getfield.(load(joinpath(@__DIR__, "helix_nebula1200px.jpg")), :r)
suite = BenchmarkGroup()

suite["helix"] = @benchmarkable ripserer(CubicalFiltration($helix))
suite["hubble"] = @benchmarkable ripserer(CubicalFiltration($hubble))

end
BenchImage.suite
