module BenchImage
using BenchmarkTools
using FileIO
using Ripserer
hubble = getfield.(load(joinpath(@__DIR__, "hubble1024px.jpg")), :r)
helix = getfield.(load(joinpath(@__DIR__, "helix_nebula1200px.jpg")), :r)
suite = BenchmarkGroup()

suite["helix"] = @benchmarkable ripserer(Cubical($helix))
suite["hubble"] = @benchmarkable ripserer(Cubical($hubble))

end
BenchImage.suite
