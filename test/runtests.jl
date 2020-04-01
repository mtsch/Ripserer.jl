using LinearAlgebra
using SparseArrays
using Test

using DataStructures

using Ripserer

function rand_dist_matrix(n)
    A = rand(n, n)
    A .+= A'
    A -= Diagonal(A)
    A
end
function rand_dist_matrix(n, sparse)
    A = sprand(n, n, sparse/2)
    A .+= A'
    A -= Diagonal(A)
    dropzeros!(A)
    A
end

@testset "Ripserer" begin
    include("helpers.jl")
    include("simplices.jl")
    include("columns.jl")
#    include("reduction.jl")
end
