"""
    copy_edges!(buffer::AbstractVector{DiameterSimplex{M}}, dist, binomial)

Copy all edges sorted by decreasing length and increasing index in distance matrix `dist` to
`buffer`.
"""
function copy_edges!(buff::AbstractVector{DiameterSimplex{M, T}},
                     dist::AbstractMatrix, binomial) where {M, T}
    empty!(buff)
    n = size(dist, 1)
    for j in 1:n, i in j+1:n
        push!(buff, DiameterSimplex{M}(T(dist[i, j]), index((i, j), binomial), 1))
    end
    sort!(buff, lt=DiameterSimplexComparer())
end

function copy_edges!(buff::AbstractVector{DiameterSimplex{M, T}},
                     dist::SparseMatrixCSC, binomial) where {M, T}
    empty!(buff)
    I, J, V = findnz(dist)
    for idx in eachindex(I)
        i = I[idx]
        j = J[idx]
        i > j || continue
        push!(buff, DiameterSimplex{M}(T(V[idx]), index((i, j), binomial), 1))
    end
    sort!(buff, lt=DiameterSimplexComparer())
end

"""
    edge_lt(e1, e2) =

Compare edges like DiameterSimplexComparer.

* by increasing diameter,
* by decreasing combinatorial index.
"""
edge_lt((e1, i1), (e2, i2)) =
    e1 < e2 || e1 == e2 && i1 > i2

"""
    edges(dist, binomial)

Get edges in distance matrix `dist`, sorted by decresing length and increasing index.
"""
function edges(dist::AbstractMatrix{T}) where T
    n = size(dist, 1)
    res = Tuple{T, Tuple{Int, Int}}[]
    for j in 1:n, i in j+1:n
        push!(res, (dist[i, j], (i, j)))
    end
    sort!(res, lt=edge_lt)
end

function edges(dist::AbstractSparseMatrix{T}) where T
    res = Tuple{T, Tuple{Int, Int}}[]
    I, J, V = findnz(dist)
    for (i, j) in zip(I, J)
        i > j || continue
        push!(res, (dist[i, j], (i, j)))
    end
    sort!(res, lt=edge_lt)
end

"""
    is_distance_matrix(dist)

Return true if dist is a valid distance matrix.
"""
is_distance_matrix(dist) =
    issymmetric(dist) && all(iszero(dist[i, i]) for i in 1:size(dist, 1))

"""
    apply_threshold(dist, thresh)

Convert matrix `dist` to sparse matrix with no entries larger than `thresh`.
"""
function apply_threshold(dist, thresh)
    n = size(dist, 1)
    for i in 1:n, j in i+1:n
        if dist[i, j] > thresh
            dist[i, j] = dist[j, i] = 0
        end
    end
    if dist isa SparseMatrixCSC
        dropzeros!(dist)
    else
        sparse(dist)
    end
end
