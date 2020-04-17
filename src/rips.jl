# distance matrix stuff ================================================================== #
"""
    edge_lt(e1, e2) =

Compare edges like DiameterSimplexComparer.

* by increasing diameter,
* by decreasing combinatorial index.
"""
edge_lt((e1, i1), (e2, i2)) =
    e1 < e2 || e1 == e2 && i1 > i2

function edges(dist::AbstractMatrix{T}, thresh=typemax(T)) where T
    n = size(dist, 1)
    res = Tuple{T, Tuple{Int, Int}}[]
    @inbounds for j in 1:n, i in j+1:n
        l = dist[i, j]
        l ≤ thresh && push!(res, (l, (i, j)))
    end
    sort!(res, lt=edge_lt)
end

function edges(dist::AbstractSparseMatrix{T}, thresh=typemax(T)) where T
    res = Tuple{T, Tuple{Int, Int}}[]
    I, J, V = findnz(dist)
    for (i, j, l) in zip(I, J, V)
        i > j || continue
        l ≤ thresh && push!(res, (l, (i, j)))
    end
    sort!(res, lt=edge_lt)
end

"""
    is_distance_matrix(dist)

Return true if dist is a valid distance matrix.
"""
is_distance_matrix(dist) =
    issymmetric(dist) && all(iszero(dist[i, i]) for i in 1:size(dist, 1))

# binomial table ========================================================================= #
"""
    Binomial(n_max, k_max)

Table of precomputed binomial coefficients up to `n_max` and `k_max`. Can be called like a
function and should be identical to [`Base.binomial`](@ref) for values of `0 ≤ n ≤ n_max`
and `0 ≤ k ≤ k_max`
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
@propagate_inbounds (bin::Binomial)(n, k) =
    bin.table[n+1, k+1]

# rips complex =========================================================================== #
"""
    default_threshold(dists)

The default threshold is equal to the radius of the input space. At this threshold, all
vertices are connected to a vertex `x` and the homology becomes trivial.
"""
default_threshold(dists) =
    minimum(maximum(dists[:, i]) for i in 1:size(dists, 1))

"""
    RipsFiltration{T, S<:AbstractSimplex{<:Any, T}}

This type holds the information about the input values.
The distance matrix has to be a dense matrix.

# Constructor

    RipsFiltration(distance_matrix;
                   dim_max=1,
                   modulus=2,
                   threshold=default_threshold(dist),
                   eltype=Simplex{modulus, T})
"""
struct RipsFiltration{T, S<:AbstractSimplex{<:Any, T}, A<:AbstractArray{T}}<:
    AbstractFiltration{T, S}

    dist         ::A
    binomial     ::Binomial
    dim_max      ::Int
    threshold    ::T
    vertex_cache ::Vector{Int}
end

function RipsFiltration(dist::AbstractArray{T};
                        dim_max::Integer=1,
                        modulus=2,
                        threshold=default_threshold(dist),
                        eltype::DataType=Simplex{modulus, T}) where T

    is_distance_matrix(dist) ||
        throw(ArgumentError("`dist` must be a distance matrix"))
    isprime(modulus) ||
        throw(ArgumentError("`modulus` must be prime"))
    dim_max ≥ 0 ||
        throw(ArgumentError("`dim_max` must be non-negative"))
    eltype <: AbstractSimplex{<:Any, T} ||
        throw(ArgumentError("`eltype` must be a subtype of `AbstractSimplex`"))
    !issparse(dist) ||
        throw(ArgumentError("`dist` is sparse. Use `SparseRipsFiltration` instead"))
    RipsFiltration{T, eltype, typeof(dist)}(
        dist, Binomial(size(dist, 1), dim_max+2), dim_max, T(threshold), Int[])
end

Base.length(rips::RipsFiltration) =
    size(rips.dist, 1)

@propagate_inbounds dist(rips::RipsFiltration, i::Integer, j::Integer) =
    rips.dist[i, j]

@propagate_inbounds Base.binomial(rips::RipsFiltration, n, k) =
    rips.binomial(n, k)

edges(rips::RipsFiltration) =
    edges(rips.dist, threshold(rips))

dim_max(rips::RipsFiltration) =
    rips.dim_max

threshold(rips::RipsFiltration) =
    rips.threshold

"""
    SparseRipsFiltration{T, S<:AbstractSimplex{<:Any, T}}

This type holds the information about the input values.
The distance matrix will be converted to a sparse matrix with all values greater than
threshold deleted. Off-diagonal zeros in the matrix are treaded as `typemax(T)`.

# Constructor

    SparseRipsFiltration(distance_matrix;
                         dim_max=1,
                         modulus=2,
                         threshold=default_threshold(dist),
                         eltype=Simplex{modulus, T})
"""
struct SparseRipsFiltration{T, S<:AbstractSimplex{<:Any, T}, A<:AbstractSparseArray{T}}<:
    AbstractFiltration{T, S}

    dist         ::A
    binomial     ::Binomial
    dim_max      ::Int
    threshold    ::T
    vertex_cache ::Vector{Int}
end

function SparseRipsFiltration(dist::AbstractArray{T};
                              dim_max::Integer=1,
                              modulus=2,
                              threshold=default_threshold(dist),
                              eltype::DataType=Simplex{modulus, T}) where T

    is_distance_matrix(dist) ||
        throw(ArgumentError("`dist` must be a distance matrix"))
    isprime(modulus) ||
        throw(ArgumentError("`modulus` must be prime"))
    dim_max ≥ 0 ||
        throw(ArgumentError("`dim_max` must be non-negative"))
    eltype <: AbstractSimplex{<:Any, T} ||
        throw(ArgumentError("`eltype` must be a subtype of `AbstractSimplex`"))

    # We need to make a copy beacuse we're editing the matrix.
    new_dist = sparse(dist)
    SparseArrays.fkeep!(new_dist, (_, _ , v) -> v ≤ threshold)

    SparseRipsFiltration{T, eltype, typeof(new_dist)}(
        new_dist, Binomial(size(dist, 1), dim_max+2), dim_max, T(threshold), Int[])
end

Base.length(rips::SparseRipsFiltration) =
    size(rips.dist, 1)

@propagate_inbounds function dist(rips::SparseRipsFiltration{T},
                                  i::Integer, j::Integer) where T
    res = rips.dist[i, j]
    ifelse(i == j, zero(T), ifelse(iszero(res), typemax(T), res))
end

@propagate_inbounds Base.binomial(rips::SparseRipsFiltration, n, k) =
    rips.binomial(n, k)

edges(rips::SparseRipsFiltration) =
    edges(rips.dist, threshold(rips))

dim_max(rips::SparseRipsFiltration) =
    rips.dim_max

threshold(rips::SparseRipsFiltration) =
    rips.threshold

SparseArrays.issparse(::Type{<:SparseRipsFiltration}) =
    true
