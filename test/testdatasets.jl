# This file includes various datasets and dataset generators used in testing.
using Distances
using LinearAlgebra
using SparseArrays
using Random

"""
    rand_dist_matrix(n, [sparse])

Construct a random distance matrix with `n` rows and columns.
"""
function rand_dist_matrix(n)
    A = rand(n, n)
    A .+= A'
    A -= Diagonal(A)
    return A
end
function rand_dist_matrix(n, sparse)
    A = sprand(n, n, sparse / 2)
    A .+= A'
    A -= Diagonal(A)
    dropzeros!(A)
    return A
end

"Distances on an icosahedron graph with edge length 1."
icosahedron = Float64[
    0 1 2 2 1 2 1 1 2 2 1 3
    1 0 3 2 1 1 2 1 2 1 2 2
    2 3 0 1 2 2 1 2 1 2 1 1
    2 2 1 0 3 2 1 1 2 1 2 1
    1 1 2 3 0 1 2 2 1 2 1 2
    2 1 2 2 1 0 3 2 1 1 2 1
    1 2 1 1 2 3 0 1 2 2 1 2
    1 1 2 1 2 2 1 0 3 1 2 2
    2 2 1 2 1 1 2 3 0 2 1 1
    2 1 2 1 2 1 2 1 2 0 3 1
    1 2 1 2 1 2 1 2 1 3 0 2
    3 2 1 1 2 1 2 2 1 1 2 0
]

"Distances on a cycle graph with edge length 1."
cycle = [
    0 1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2 1
    1 0 1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2
    2 1 0 1 2 3 4 5 6 7 8 9 8 7 6 5 4 3
    3 2 1 0 1 2 3 4 5 6 7 8 9 8 7 6 5 4
    4 3 2 1 0 1 2 3 4 5 6 7 8 9 8 7 6 5
    5 4 3 2 1 0 1 2 3 4 5 6 7 8 9 8 7 6
    6 5 4 3 2 1 0 1 2 3 4 5 6 7 8 9 8 7
    7 6 5 4 3 2 1 0 1 2 3 4 5 6 7 8 9 8
    8 7 6 5 4 3 2 1 0 1 2 3 4 5 6 7 8 9
    9 8 7 6 5 4 3 2 1 0 1 2 3 4 5 6 7 8
    8 9 8 7 6 5 4 3 2 1 0 1 2 3 4 5 6 7
    7 8 9 8 7 6 5 4 3 2 1 0 1 2 3 4 5 6
    6 7 8 9 8 7 6 5 4 3 2 1 0 1 2 3 4 5
    5 6 7 8 9 8 7 6 5 4 3 2 1 0 1 2 3 4
    4 5 6 7 8 9 8 7 6 5 4 3 2 1 0 1 2 3
    3 4 5 6 7 8 9 8 7 6 5 4 3 2 1 0 1 2
    2 3 4 5 6 7 8 9 8 7 6 5 4 3 2 1 0 1
    1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2 1 0
]

"Projective plane, taken from ripser/examples."
projective_plane = [
    0 1 1 1 1 1 1 1 1 2 2 2 2
    1 0 2 2 2 1 2 1 2 1 2 2 2
    1 2 0 2 2 2 1 2 1 1 2 2 2
    1 2 2 0 2 1 2 2 1 2 2 2 1
    1 2 2 2 0 2 1 1 2 2 2 2 1
    1 1 2 1 2 0 2 2 2 1 1 2 1
    1 2 1 2 1 2 0 2 2 1 1 2 1
    1 1 2 2 1 2 2 0 2 1 2 1 1
    1 2 1 1 2 2 2 2 0 1 2 1 1
    2 1 1 2 2 1 1 1 1 0 1 1 2
    2 2 2 2 2 1 1 2 2 1 0 2 1
    2 2 2 2 2 2 2 1 1 1 2 0 1
    2 2 2 1 1 1 1 1 1 2 1 1 0
]

"""
    rand_n_sphere(n, dim)

Randomly sample `n` points from `dim`-dimensional sphere. Return euclidean distance matrix.
"""
function rand_n_sphere(n, dim)
    points = mapslices(normalize, randn(dim + 1, n); dims=1)
    return pairwise(Euclidean(), points)
end

"""
    torus_dist(pts)

Calculate euclidean distances between points on [-1,1]×[-1,1] as if they were on a flat
torus.
"""
function torus_dist(pts)
    n = size(pts, 2)
    dist = fill(Inf, (n, n))
    for i in 1:n
        for j in (i + 1):n
            for x_offset in -2:2:2, y_offset in -2:2:2
                off = [x_offset, y_offset]
                dst = evaluate(Euclidean(), pts[:, i], pts[:, j] .+ off)
                dist[i, j] = dist[j, i] = min(dist[i, j], dst)
            end
        end
        dist[i, i] = 0
    end
    return dist
end

"""
    torus_points(n, R=2, r=1)

Get `n` points from torus with external radius `R` and internal radius `r` embedded in 3d
euclidean space.
"""
function torus_points(n, R=2, r=1)
    n = floor(Int, sqrt(n))
    return [
        ((R + r * cos(θ)) * cos(φ), (R + r * cos(θ)) * sin(φ), r * sin(θ)) for
        θ in range(0, 2π; length=n + 1)[1:(end - 1)] for
        φ in range(0, 2π; length=n + 1)[1:(end - 1)]
    ]
end

"""
    torus(n)

Construct a torus distance matrix with `n` points. The points are equidistant. `n` must be a
square number.
"""
function torus(n)
    sqn = √n
    @assert sqn == floor(Int, sqn)
    n = floor(Int, sqn)
    r = range(-1, 1; length=n + 1)[1:(end - 1)]
    pts = fill(0.0, (2, n * n))
    i = 1
    for x in r, y in r
        pts[1, i] = x
        pts[2, i] = y
        i += 1
    end
    return torus_dist(pts)
end

"""
    rand_torus(n)

Construct a random torus distance matrix with `n` points.
"""
function rand_torus(n)
    pts = rand(2, n) .* 2 .- 1
    return torus_dist(pts)
end

"""
    disconnected_tori(n, m)

Construct sparse distance matrix of `m` random toruses with `n` points each.
"""
function disconnected_tori(n, m)
    dist = zeros(n * m, n * m)
    for i in 1:n:(n * m)
        dist[i:(i + n - 1), i:(i + n - 1)] .= rand_torus(n)
    end
    perm = shuffle(1:(n * m))
    return sparse(dist[perm, perm])
end
