using LinearAlgebra
using SparseArrays
using Test

using DataStructures

using Ripserer

@testset "Ripserer" begin
    include("helpers.jl")
    include("simplices.jl")
    include("distancematrix.jl")
    include("columns.jl")
    include("reduction.jl")
end
