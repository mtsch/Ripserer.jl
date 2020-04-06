using LinearAlgebra
using SparseArrays
using Test

using DataStructures

using Ripserer
include("data.jl")

@testset "Ripserer" begin
    include("rips.jl")
    include("reduction.jl")
end
