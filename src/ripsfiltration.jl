# distance matrix stuff ================================================================== #
"""
    edges(dist::AbstractMatrix{T}, thresh, E)

Return sorted edges of type `E` in distance matrix with length lower than `thresh`.
"""
function edges(dist::AbstractMatrix{T}, thresh, edge_type) where T
    @boundscheck begin
        size(dist, 1) == size(dist, 2) && size(dist, 1) > 1 ||
            throw(ArgumentError("invalid input matrix"))
    end

    n = size(dist, 1)
    # infer type
    E = typeof(edge_type((2, 1), dist[2, 1]))
    res = E[]
    @inbounds for j in 1:n, i in j+1:n
        l = dist[i, j]
        l ≤ thresh && push!(res, E((i, j), l))
    end
    sort!(res)
end

function edges(dist::AbstractSparseMatrix{T}, thresh, edge_type) where T
    # infer type
    E = typeof(edge_type((2, 1), dist[2, 1]))
    res = E[]
    I, J, V = findnz(dist)
    for (i, j, l) in zip(I, J, V)
        i > j || continue
        l ≤ thresh && push!(res, E((i, j), l))
    end
    sort!(res)
end

"""
    distances(metric, points[, births])

Return distance matrix calculated from `points` with `metric`. If given, add birth times
from `births` to the diagonal.
"""
function distances(metric, points, births=nothing)
    dim = length(first(points))
    T = eltype(first(points))
    dists = pairwise(metric, reshape(reinterpret(T, points), (dim, length(points))), dims=2)
    if !isnothing(births)
        for i in 1:length(points)
            dists[i, i] = births[i]
        end
    end
    dists
end

# flag filtration ======================================================================== #
"""
    AbstractFlagFiltration{T, V} <: AbstractFiltration{T, V}

An abstract flag filtration is a filtration of flag complexes. Its subtypes can overload
`dist(::AbstractFlagFiltration{T}, u, v)::Union{T, Infinity}` instead of `diam`.
`diam(::AbstractFlagFiltration, ...)` defaults to maximum `dist` among vertices.
"""
abstract type AbstractFlagFiltration{T, V} <: AbstractFiltration{T, V} end

@propagate_inbounds function diam(flt::AbstractFlagFiltration, vertices)
    n = length(vertices)
    res = typemin(dist_type(flt))
    for i in 1:n, j in i+1:n
        d = dist(flt, vertices[j], vertices[i])
        res = ifelse(res > d, res, d)
    end
    ifelse(res > threshold(flt), ∞, res)
end

@propagate_inbounds function diam(flt::AbstractFlagFiltration, sx::AbstractSimplex, us, v)
    res = diam(sx)
    for u in us
        # Even though this looks like a tight loop, v changes way more often than us, so
        # this is the faster order of indexing by u and v.
        d = dist(flt, v, u)
        res = ifelse(res > d, res, d)
    end
    ifelse(res > threshold(flt), ∞, res)
end

edges(flt::AbstractFlagFiltration) =
    edges(flt.dist, threshold(flt), edge_type(flt))

"""
    dist(::AbstractFlagFiltration, u, v)

Return the distance between vertices `u` and `v`. If the distance is higher than the
threshold, return `Infinity()` instead.
"""
dist(::AbstractFlagFiltration, ::Any, ::Any)

# rips filtration ======================================================================== #
"""
    default_rips_threshold(dists)

The default threshold is equal to the radius of the input space. At this threshold, all
vertices are connected to a vertex `x` and the homology becomes trivial. If any of the
distances is negative, default threshold defaults to `typemax(eltype(dists))`.
"""
default_rips_threshold(dists::AbstractMatrix{T}) where T =
    minimum(maximum(abs, dists[:, i]) for i in 1:size(dists, 1))

"""
    RipsFiltration{T, V<:AbstractSimplex{<:Any, T}} <: AbstractFlagFiltration{T, V}

This type represents a filtration of Vietoris-Rips complexes.
Diagonal items are treated as vertex birth times.

# Constructor

    RipsFiltration(
        distance_matrix;
        threshold=default_rips_threshold(dist),
        vertex_type=Simplex{0, T, Int64},
    )
"""
struct RipsFiltration{
    T, V<:AbstractSimplex{0, T}, A<:AbstractMatrix{T}
} <: AbstractFlagFiltration{T, V}

    dist      ::A
    threshold ::T
end

function RipsFiltration(
    dist::AbstractMatrix{T};
    threshold=default_rips_threshold(dist),
    vertex_type::DataType=Simplex{0, T, Int64}
) where T

    issymmetric(dist) ||
        throw(ArgumentError("`dist` must be a distance matrix"))
    vertex_type <: AbstractSimplex{0, T} ||
        throw(ArgumentError("`vertex_type` must be a subtype of `AbstractSimplex{0, $T}`"))
    !issparse(dist) ||
        throw(ArgumentError("`dist` is sparse. Use `SparseRipsFiltration` instead"))
    RipsFiltration{T, vertex_type, typeof(dist)}(dist, T(threshold))
end

# interface implementation
n_vertices(rips::RipsFiltration) =
    size(rips.dist, 1)

@propagate_inbounds dist(rips::RipsFiltration{T}, i, j) where T =
    ifelse(i == j, zero(T), rips.dist[i, j])

threshold(rips::RipsFiltration) =
    rips.threshold

birth(rips::RipsFiltration, i) =
    rips.dist[i, i]

# sparse rips filtration ================================================================= #
"""
    SparseRipsFiltration{T, V<:AbstractSimplex{<:Any, T}} <: AbstractFlagFiltration{T, V}

This type represents a filtration of Vietoris-Rips complexes.
The distance matrix will be converted to a sparse matrix with all values greater than
threshold deleted. Off-diagonal zeros in the matrix are treated as ∞. Diagonal items are
treated as vertex birth times.

# Constructor

    SparseRipsFiltration(
        distance_matrix;
        threshold=nothing,
        vertex_type=Simplex{0, T, Int64},
    )
"""
struct SparseRipsFiltration{
    T, V<:AbstractSimplex{0, T}, A<:AbstractSparseMatrix{T}
}<: AbstractFlagFiltration{T, V}

    dist ::A
end

function SparseRipsFiltration(
    dist::AbstractMatrix{T};
    threshold=nothing,
    vertex_type::DataType=Simplex{0, T, Int64},
) where T
    issymmetric(dist) ||
        throw(ArgumentError("`dist` must be a distance matrix"))
    vertex_type <: AbstractSimplex{0, T} ||
        throw(ArgumentError("`vertex_type` must be a subtype of `AbstractSimplex{0, $T}`"))

    new_dist = sparse(dist)
    if !isnothing(threshold)
        births = diag(dist)
        SparseArrays.fkeep!(new_dist, (_, _ , v) -> v ≤ threshold)
        for i in 1:size(dist, 1)
            new_dist[i, i] = births[i]
        end
    end
    SparseRipsFiltration{T, vertex_type, typeof(new_dist)}(new_dist)
end

# interface implementation --------------------------------------------------------------- #
n_vertices(rips::SparseRipsFiltration) =
    size(rips.dist, 1)

@propagate_inbounds function dist(rips::SparseRipsFiltration{T}, i, j) where T
    res = rips.dist[i, j]
    ifelse(i == j, zero(T), ifelse(iszero(res), ∞, res))
end

# Threshold was handled by deleting entries in the matrix.
threshold(rips::SparseRipsFiltration) =
    ∞

@propagate_inbounds function diam(rips::SparseRipsFiltration, sx::AbstractSimplex, us, v)
    res = diam(sx)
    for u in us
        # Since indexing in sparse matrices is expensive, we want to abort the loop early
        # even though the number of vertices in us is small.
        d = dist(rips, v, u)
        d == ∞ && return ∞
        res = ifelse(res > d, res, d)
    end
    res
end

birth(rips::SparseRipsFiltration, i) =
    rips.dist[i, i]
