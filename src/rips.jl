# basic simplex ========================================================================== #
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
    n_bits(M)

Get numer of bits needed to represent number mod `M`.
"""
n_bits(M) =
    floor(Int, log2(M-1)) + 1

"""
    Simplex{M, T} <: AbstractSimplex{M, T}

The basic simplex type with coefficient values from `Z_M`, integers modulo `M`.
`index` and `coef` are packed into a single `UInt64`.

# Constructor

    Simplex{M}(::T, index::Integer, coef::Integer)
"""
struct Simplex{M, T} <: AbstractSimplex{M, T}
    diam       ::T
    index_coef ::UInt64

    # Prevent accidentally creating simplex with Simplex{M, T}(::T, ::Int)
    Simplex{M, T}(diam::T, index_coef::UInt64) where {M, T} =
        new{M, T}(diam, index_coef)
end

@generated function Simplex{M}(diam::T, index, coef) where {M, T}
    isprime(M) || throw(DomainError(M, "modulus not prime"))
    bits = n_bits(M)
    :(Simplex{M, T}(diam, UInt64(index) << $bits + mod(coef, $M)))
end
Simplex{M, T}(diam::T, index, coef) where {M, T} =
    Simplex{M}(diam, index, coef)
Simplex{M}(scx::SimplicialComplex{M}, diam, vertices, coef) where M =
    Simplex{M}(diam, index(scx, vertices), coef)

@generated function index(sx::Simplex{M}) where M
    shift = n_bits(M)
    :(reinterpret(Int64, sx.index_coef >> $shift))
end

@generated function coef(sx::Simplex{M}) where M
    mask = 1 << n_bits(M) - 1
    :(reinterpret(Int64, sx.index_coef & $mask))
end

diam(sx::Simplex) =
    sx.diam

@generated function set_coef(sx::Simplex{M, T}, value) where {M, T}
    mask = ~(1 << n_bits(M) - 1)
    :(Simplex{M, T}(diam(sx), sx.index_coef & $mask + mod(value, M)))
end

Base.show(io::IO, sx::Simplex{M}) where M =
    print(io, "Simplex{", M, "}", (diam(sx), index(sx), coef(sx)))

# distance matrix stuff ================================================================== #
"""
    edge_lt(e1, e2) =

Compare edges like DiameterSimplexComparer.

* by increasing diameter,
* by decreasing combinatorial index.
"""
edge_lt((e1, i1), (e2, i2)) =
    e1 < e2 || e1 == e2 && i1 > i2

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

# binomial table ========================================================================= #
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

# rips complex =========================================================================== #
"""
    RipsComplex{M, T, A<:AbstractArray{T}}

This type holds the information about the input values.

# Constructor

    RipsComplex{M}(distance_matrix, dim_max)

"""
struct RipsComplex{M, T, A<:AbstractArray{T}} <: SimplicialComplex{M, T, Simplex{M, T}}
    dist         ::A
    binomial     ::Binomial
    dim_max      ::Int
    vertex_cache ::Vector{Int}

    function RipsComplex{M}(dist::A, dim_max::Integer) where {M, T, A<:AbstractArray{T}}
        is_distance_matrix(dist) ||
            throw(ArgumentError("`dist` must be a distance matrix"))
        isprime(M) ||
            throw(ArgumentError("`modulus` must be prime"))
        dim_max ≥ 0 ||
            throw(ArgumentError("`dim_max` must be positive"))
        new{M, T, A}(dist, Binomial(size(dist, 1), dim_max+2), dim_max, Int[])
    end
end

Base.length(rips::RipsComplex) =
    size(rips.dist, 1)

dist(rips::RipsComplex, i::Integer, j::Integer) =
    rips.dist[i, j]

Base.binomial(rips::RipsComplex, n, k) =
    rips.binomial(n, k)

edges(rips::RipsComplex) =
    edges(rips.dist)

dim_max(rips::RipsComplex) =
    rips.dim_max
