using Compat
using LinearAlgebra
using SparseArrays
using Test

using DataStructures

using Ripserer
include("data.jl")

using Aqua

@testset "Ripserer" begin
    Aqua.test_all(Ripserer)

    include("simplex.jl")
    include("rips.jl")
    include("reduction.jl")
end
