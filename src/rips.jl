"""
    AbstractFlagFiltration{T, S} <: AbstractFiltration{T, S}

An abstract flag filtration is a filtration of flag complexes. Its subtypes can overload
`dist(::AbstractFlagFiltration{T}, u, v)::Union{T, Infinity}` instead of `diam`.
`diam(::AbstractFlagFiltration, ...)` defaults to maximum `dist` among vertices.
"""
abstract type AbstractFlagFiltration{T, S} <: AbstractFiltration{T, S} end

@propagate_inbounds function diam(flt::AbstractFlagFiltration, vertices)
    n = length(vertices)
    res = typemin(dist_type(flt))
    for i in 1:n, j in i+1:n
        d = dist(flt, vertices[j], vertices[i])
        res = ifelse(res > d, res, d)
    end
    if res > threshold(flt)
        ∞
    else
        res
    end
end

edges(flt::AbstractFlagFiltration) =
    edges(flt.dist, threshold(flt))

"""
    dist(::AbstractFlagFiltration, u, v)

Return the distance between vertices `u` and `v`. If the distance is higher than the
threshold, return `Infinity()` instead.
"""
dist

"""
    default_rips_threshold(dists)

The default threshold is equal to the radius of the input space. At this threshold, all
vertices are connected to a vertex `x` and the homology becomes trivial.
"""
default_rips_threshold(dists) =
    minimum(maximum(dists[:, i]) for i in 1:size(dists, 1))

"""
    RipsFiltration{T, S<:AbstractSimplex{<:Any, T}} <: AbstractFlagFiltration{T, S}

# Constructor

    RipsFiltration(distance_matrix;
                   modulus=2,
                   threshold=default_rips_threshold(dist),
                   simplex_type=Simplex{modulus, T})
"""
struct RipsFiltration{
    T, S<:AbstractSimplex{<:Any, T}, A<:AbstractMatrix{T}
} <: AbstractFlagFiltration{T, S}

    dist      ::A
    threshold ::T
end

function RipsFiltration(
    dist::AbstractMatrix{T};
    modulus=2,
    threshold=default_rips_threshold(dist),
    simplex_type::DataType=Simplex{modulus, T}
) where T

    is_distance_matrix(dist) ||
        throw(ArgumentError("`dist` must be a distance matrix"))
    is_prime(modulus) ||
        throw(ArgumentError("`modulus` must be prime"))
    simplex_type <: AbstractSimplex{<:Any, T} ||
        throw(ArgumentError("`simplex_type` must be a subtype of `AbstractSimplex`"))
    !issparse(dist) ||
        throw(ArgumentError("`dits` is sparse. Use `SparseRipsFiltration` instead"))

    if !isfinite(threshold)
        threshold = default_rips_threshold(dist)
    end
    RipsFiltration{T, simplex_type, typeof(dist)}(dist, T(threshold))
end
RipsFiltration(points; metric=Euclidean(), kwargs...) =
    RipsFiltration(distances(metric, points); kwargs...)

n_vertices(rips::RipsFiltration) =
    size(rips.dist, 1)

@propagate_inbounds dist(rips::RipsFiltration, i, j) =
    rips.dist[i, j]

threshold(rips::RipsFiltration) =
    rips.threshold

@propagate_inbounds function diam(flt::AbstractFlagFiltration, sx, us, v::Integer)
    res = diam(sx)
    for u in us
        # Even though this looks like a tight loop, v changes way more often than us, so
        # this is the faster order of indexing by u and v.
        d = dist(flt, v, u)
        res = ifelse(res > d, res, d)
    end
    if res > threshold(flt)
        ∞
    else
        res
    end
end

"""
    SparseRipsFiltration{T, S<:AbstractSimplex{<:Any, T}} <: AbstractFlagFiltration{T, S}

This type holds the information about the input values.
The distance matrix will be converted to a sparse matrix with all values greater than
threshold deleted. Off-diagonal zeros in the matrix are treated as `typemax(T)`.

# Constructor

    SparseRipsFiltration(distance_matrix;
                         modulus=2,
                         threshold=default_rips_threshold(dist),
                         eltype=Simplex{modulus, T})
"""
struct SparseRipsFiltration{
    T, S<:AbstractSimplex{<:Any, T}, A<:AbstractSparseMatrix{T}
}<: AbstractFlagFiltration{T, S}

    dist         ::A
    threshold    ::T
end

function SparseRipsFiltration(dist::AbstractMatrix{T};
                              modulus=2,
                              threshold=default_rips_threshold(dist),
                              simplex_type::DataType=Simplex{modulus, T}) where T

    is_distance_matrix(dist) ||
        throw(ArgumentError("`dist` must be a distance matrix"))
    is_prime(modulus) ||
        throw(ArgumentError("`modulus` must be prime"))
    simplex_type <: AbstractSimplex{<:Any, T} ||
        throw(ArgumentError("`simplex_type` must be a subtype of `AbstractSimplex`"))
    if !isfinite(threshold)
        threshold = default_rips_threshold(dist)
    end
    # We need to make a copy beacuse we're editing the matrix.
    new_dist = sparse(dist)
    SparseArrays.fkeep!(new_dist, (_, _ , v) -> v ≤ threshold)

    SparseRipsFiltration{T, simplex_type, typeof(new_dist)}(new_dist, T(threshold))
end
SparseRipsFiltration(points; metric=Euclidean(), kwargs...) =
    SparseRipsFiltration(distances(metric, points); kwargs...)

n_vertices(rips::SparseRipsFiltration) =
    size(rips.dist, 1)

@propagate_inbounds function dist(rips::SparseRipsFiltration{T}, i, j) where T
    res = rips.dist[i, j]
    ifelse(i == j, zero(T), ifelse(iszero(res), ∞, res))
end

@propagate_inbounds Base.binomial(rips::SparseRipsFiltration, n, k) =
    rips.binomial(n, k)

# Threshold was handled by deleting entries in the matrix.
threshold(rips::SparseRipsFiltration) =
    ∞

SparseArrays.issparse(::Type{<:SparseRipsFiltration}) =
    true

@propagate_inbounds function diam(flt::SparseRipsFiltration, sx, us, v::Integer)
    res = diam(sx)
    for u in us
        # Since indexing in sparse matrices is expensive, we want to abort the loop early
        # even though the number of vertices in us is small.
        d = dist(flt, v, u)
        d == ∞ && return ∞
        res = ifelse(res > d, res, d)
    end
    res
end
