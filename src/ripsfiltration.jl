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
    return sort!(res)
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
    return sort!(res)
end

"""
    distances(metric, points[, births])

Return distance matrix calculated from `points` with `metric`. If given, add birth times
from `births` to the diagonal.
"""
function distances(metric, points)
    isempty(points) && throw(ArgumentError("`points` must be nonempty"))

    dim = length(first(points))
    T = eltype(first(points))
    dists = pairwise(metric, reshape(reinterpret(T, points), (dim, length(points))), dims=2)
    return dists
end

# rips filtrations ======================================================================= #
"""
    AbstractRipsFiltration{T, S} <: AbstractFiltration{T, S}

An abstract Vietoris-Rips filtration. Its subtypes can overload
`dist(::AbstractRipsFiltration{T}, u, v)::Union{T, Missing}` instead of `diam`.
`diam(::AbstractRipsFiltration, ...)` defaults to maximum `dist` among vertices.
"""
abstract type AbstractRipsFiltration{T, S} <: AbstractFiltration{T, S} end

@propagate_inbounds function diam(flt::AbstractRipsFiltration, vertices)
    n = length(vertices)
    res = typemin(dist_type(flt))
    for i in 1:n, j in i+1:n
        d = dist(flt, vertices[j], vertices[i])
        ismissing(d) && return missing
        res = ifelse(d < res, res, d)
    end
    return ifelse(res > threshold(flt), missing, res)
end

@propagate_inbounds function diam(flt::AbstractRipsFiltration, sx::AbstractSimplex, us, v)
    res = diam(sx)
    for u in us
        # Even though this looks like a tight loop, v changes way more often than us, so
        # this is the faster order of indexing by u and v.
        d = dist(flt, v, u)
        if ismissing(d) || d > threshold(flt)
            return missing
        else
            res = ifelse(res > d, res, d)
        end
    end
    return res
end

edges(flt::AbstractRipsFiltration) = edges(flt.dist, threshold(flt), edge_type(flt))

"""
    dist(::AbstractRipsFiltration, u, v)

Return the distance between vertices `u` and `v`. If the distance is higher than the
threshold, return `missing` instead.
"""
dist(::AbstractRipsFiltration, ::Any, ::Any)

# rips filtration ======================================================================== #
"""
    default_rips_threshold(dists)

The default threshold is equal to the radius of the input space. At this threshold, all
vertices are connected to a vertex `x` and the homology becomes trivial. If any of the
distances is negative, default threshold defaults to `typemax(eltype(dists))`.
"""
function default_rips_threshold(dists::AbstractMatrix{T}) where T
    return minimum(maximum(abs, dists[:, i]) for i in 1:size(dists, 1))
end

"""
    Rips{I, T} <: AbstractRipsFiltration{T, Simplex}

This type represents a filtration of Vietoris-Rips complexes.
Diagonal items are treated as vertex birth times.

# Constructors

* `Rips(distance_matrix; threshold=default_rips_threshold(dist))`
* `Rips(points; metric=Euclidean(), threshold=default_rips_threshold(dist))`
* `Rips{I}(args...)`: `I` sets the size of integer used to represent simplices.
"""
struct Rips{I<:Integer, T, A<:AbstractMatrix{T}} <: AbstractRipsFiltration{T, Simplex}
    dist::A
    threshold::T
end

function Rips{I}(
    dist::AbstractMatrix{T}; threshold=nothing
) where {I, T}
    issymmetric(dist) ||
        throw(ArgumentError("`dist` must be a distance matrix"))
    !issparse(dist) ||
        throw(ArgumentError("`dist` is sparse. Use `SparseRips` instead"))
    thresh = isnothing(threshold) ? default_rips_threshold(dist) : T(threshold)
    return Rips{I, T, typeof(dist)}(dist, thresh)
end
function Rips{I}(points::AbstractVector; metric=Euclidean(), kwargs...) where I
    return Rips(distances(metric, points); kwargs...)
end
function Rips(dist; kwargs...)
    return Rips{Int}(dist; kwargs...)
end

@propagate_inbounds function dist(rips::Rips{<:Any, T}, i, j) where T
    return ifelse(i == j, zero(T), rips.dist[i, j])
end

n_vertices(rips::Rips) = size(rips.dist, 1)
threshold(rips::Rips) = rips.threshold
birth(rips::Rips, i) = rips.dist[i, i]
simplex_type(rips::Rips{I, T}, dim) where {I, T} = Simplex{dim, T, I}

# sparse rips filtration ================================================================= #
"""
    SparseRips{I, T} <: AbstractRipsFiltration{T, Simplex}

This type represents a filtration of Vietoris-Rips complexes.
The distance matrix will be converted to a sparse matrix with all values greater than
threshold deleted. Off-diagonal zeros in the matrix are treated as `missing`. Diagonal items
are treated as vertex birth times.

# Constructor

* `SparseRips{I}(distance_matrix; threshold=nothing)`
* `SparseRips(distance_matrix; threshold=nothing)`: `I` sets the integer size used to
  represent simplices.
"""
struct SparseRips{I, T, A<:AbstractSparseMatrix{T}} <: AbstractRipsFiltration{T, Simplex}
    dist::A
    threshold::T
end

function SparseRips{I}(
    dist::AbstractSparseMatrix{T}; threshold=nothing
) where {I, T}
    issymmetric(dist) || throw(ArgumentError("`dist` must be a distance matrix"))
    if isnothing(threshold)
        thresh = maximum(dist)
        dists = sparse(dist)
    else
        thresh = threshold
        dists = SparseArrays.fkeep!(copy(dist), (_, _, v) -> v ≤ threshold)
    end
    return SparseRips{I, T, typeof(dists)}(dists, thresh)
end
function SparseRips(dist; threshold=nothing)
    return SparseRips{Int}(dist; threshold=threshold)
end
function SparseRips{I}(dist; threshold=nothing) where I
    return SparseRips{I}(sparse(dist); threshold=threshold)
end

@propagate_inbounds function dist(rips::SparseRips{<:Any, T}, i, j) where T
    res = rips.dist[i, j]
    return ifelse(i == j, zero(T), ifelse(iszero(res), missing, res))
end

n_vertices(rips::SparseRips) = size(rips.dist, 1)
threshold(rips::SparseRips) = rips.threshold
birth(rips::SparseRips, i) = rips.dist[i, i]
simplex_type(rips::SparseRips{I, T}, dim) where {I, T} = Simplex{dim, T, I}
