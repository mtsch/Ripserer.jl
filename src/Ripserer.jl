module Ripserer

using DataStructures
using SparseArrays
using LinearAlgebra

include("helpers.jl")
include("simplices.jl")
include("columns.jl")
#include("reduction.jl")

function ripserer(dist, dim_max = 1)
    res = Vector{Tuple{Int, Int}}[]
    dist = rand_dist_matrix(100)
    st = ReductionState(dist, 1, 2)
    columns = DiameterSimplex{2, Float64}[]
    res[1] = compute_0_dim_pairs!(st, columns)

    for dim in 1:dim_max
        rmx = ReductionMatrices(st)
        res[2] = compute_pairs!(rmx, columns, 1)
    end

    res
end

export ripserer
end
