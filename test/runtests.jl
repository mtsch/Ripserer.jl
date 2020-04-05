using LinearAlgebra
using SparseArrays
using Test

using DataStructures

using Ripserer

"""
    rand_dist_matrix(n, [sparse])

Construct a random distance matrix with `n` rows and columns.
"""
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

"""
    torus(n)

Construct a torus distance matrix with `n` points. The points are equidistant.
"""
function torus(n)
    n = floor(Int, sqrt(n))
    r = range(-1, 1, length = n+1)[1:end-1]
    pts = [(i, j) for i in r for j in r]
    dist = fill(Inf, length(pts), length(pts))
    for i in 1:length(pts), j in i+1:length(pts)
        for x_offset in -2:2:2, y_offset in -2:2:2
            x1, y1 = pts[i]
            x2, y2 = pts[j] .+ (x_offset, y_offset)
            dist[i, j] = dist[j, i] = min(√((x1-x2)^2 + (y1-y2)^2), dist[i, j])
        end
    end
    for i in 1:length(pts)
        dist[i, i] = 0
    end
    dist
end

@testset "Ripserer" begin
    include("helpers.jl")
    include("simplices.jl")
    include("columns.jl")
end
