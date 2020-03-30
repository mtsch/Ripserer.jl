"""
    isprime(n)

Return `true` if `n` is a prime number.
"""
function isprime(n)
    if iseven(n) || n < 2
        n == 2
    else
        p = 3
        q = n / p
        while p ≤ q
            iszero(n % p) && return false
            p += 2
            q = n / p
        end
        true
    end
end

"""
    Binomial(n_max, k_max)

Table of precomputed binomial coefficients up to `n_max` and `k_max`. Can be called like a
function and should be identical to [`binomial`](@ref) for values of `0 ≤ n ≤ n_max` and
`0 ≤ k ≤ k_max`
"""
struct Binomial
    table::Matrix{Int64}
end

function Binomial(n, k)
    table = zeros(Int, n+1, k+1)
    for i in 1:n+1
	table[i, 1] = 1;
	for j in 2:min(i, k+1)
	    table[i, j] = table[i-1, j-1] + table[i-1, j];
	    if (i <= k)
                table[i, i] = 1
            end
        end
    end
    Binomial(table)
end

Base.show(io::IO, bin::Binomial) =
    print(io, "Binomial$(size(bin.table) .- 1)")
(bin::Binomial)(n, k) =
    bin.table[n+1, k+1]

# distance matrix stuff ================================================================== #
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

# compressed sparse matrix =============================================================== #
"""
    CompressedSparseMatrix{T}

Compressed immutable sparse matrix data structure that supports efficient column insertion,
pushing to the last column via [`push!`](@ref) and iterating over columns.
It's up to the value type `T` to know about its row position.
"""
struct CompressedSparseMatrix{T}
    colptr::Vector{Int}
    nzval::Vector{T}
end

CompressedSparseMatrix{T}() where T =
    CompressedSparseMatrix(Int[1], T[])

function Base.show(io::IO, csm::CompressedSparseMatrix{T}) where T
    println(io, "CompressedSparseMatrix{$T}[")
    for i in 1:length(csm)
        println(io, "  $i: ", collect(csm[i]))
    end
    print(io, "]")
end

function Base.push!(csm::CompressedSparseMatrix, value)
    push!(csm.nzval, value)
    csm.colptr[end] += 1
    value
end

add_column!(csm::CompressedSparseMatrix) =
    push!(csm.colptr, csm.colptr[end])
Base.eltype(csm::CompressedSparseMatrix{T}) where T =
    T
Base.length(csm::CompressedSparseMatrix) =
    length(csm.colptr) - 1
Base.getindex(csm::CompressedSparseMatrix, i) =
    CSMColumnIterator(csm, i)

struct CSMColumnIterator{T}
    csm ::CompressedSparseMatrix{T}
    idx ::Int
end

Base.IteratorSize(::Type{CSMColumnIterator}) =
    Base.HasLength()
Base.IteratorEltype(::Type{CSMColumnIterator{T}}) where T =
    Base.HasEltype()
Base.eltype(::Type{CSMColumnIterator{T}}) where T =
    T
Base.length(ci::CSMColumnIterator) =
    ci.csm.colptr[ci.idx + 1] - ci.csm.colptr[ci.idx]

function Base.iterate(ci::CSMColumnIterator, i=1)
    colptr = ci.csm.colptr
    index = i + colptr[ci.idx] - 1
    if index ≥ colptr[ci.idx + 1]
        nothing
    else
        (ci.csm.nzval[index], i + 1)
    end
end
