using LinearAlgebra
using SparseArrays
using Test

using DataStructures

using Ripserer
include("data.jl")

@testset "Ripserer" begin
    include("helpers.jl")
    include("simplices.jl")
    include("columns.jl")
end
