module BenchImage
using BenchmarkTools
using FileIO
using Ripserer
hubble1024 = getfield.(load(joinpath(@__DIR__, "hubble1024px.jpg")), :r)

bonsai256 = Array{UInt8}(undef, (256, 256, 256))
read!(joinpath(@__DIR__, "bonsai.raw"), bonsai256)
bonsai64 = bonsai256[1:4:end, 1:4:end, 1:4:end]

suite = BenchmarkGroup()

suite["hubble"] = @benchmarkable ripserer(Cubical($hubble1024))
suite["bonsai"] = @benchmarkable ripserer(Cubical($bonsai64))

end
BenchImage.suite
