using Compat
using LinearAlgebra
using SparseArrays
using Test

using Aqua
using DataStructures
using Distances

using Ripserer

include("data.jl")

@testset "Ripserer" begin
    Aqua.test_all(Ripserer)

    include("convenience.jl")
    include("simplex.jl")
    include("rips.jl")
    include("reduction.jl")
end
