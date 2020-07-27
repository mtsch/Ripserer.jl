module BenchImage
using BenchmarkTools
using FileIO
using Ripserer
hubble1024 = getfield.(load(joinpath(@__DIR__, "hubble1024px.jpg")), :r)

# TODO: add bonsai when supported by new version of Cubical
# bonsai256 = Array{UInt8}(undef, (256, 256, 256))
# read!(joinpath(@__DIR__, "bonsai.raw"), bonsai256)
# bonsai64 = bonsai256[1:4:end, 1:4:end, 1:4:end]

suite = BenchmarkGroup()

suite["hubble"] = @benchmarkable ripserer(Cubical{Int128}($hubble1024))

end
BenchImage.suite
