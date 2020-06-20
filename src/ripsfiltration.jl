# distance matrix stuff ================================================================== #
"""
    edges(dist::AbstractMatrix, thresh, E)

Return sorted edges of type `E` in distance matrix with length lower than `thresh`.
"""
function edges(dist::AbstractMatrix, thresh, ::Type{E}) where E
    @boundscheck begin
        size(dist, 1) == size(dist, 2) && size(dist, 1) > 1 ||
            throw(ArgumentError("invalid input matrix"))
    end

    n = size(dist, 1)
    res = E[]
    @inbounds for j in 1:n, i in j+1:n
        l = dist[i, j]
        l ≤ thresh && push!(res, E((i, j), l))
    end
    return sort!(res)
end

function edges(dist::AbstractSparseMatrix, thresh, ::Type{E}) where E
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
    AbstractRipsFiltration{I<:Signed, T} <: AbstractFiltration

An abstract Vietoris-Rips filtration. Its subtypes can overload
[`dist`](@ref)`(::AbstractRipsFiltration{T}, u, v)::Union{T, Missing}` instead of
[`diam`](@ref).

[`diam`](@ref)`(::AbstractRipsFiltration, ...)` then defaults to maximum [`dist`] among
vertices.

Comes with default implementations of [`edges`](@ref) and [`simplex_type`](@ref).
"""
abstract type AbstractRipsFiltration{I<:Signed, T} <: AbstractFiltration end

@propagate_inbounds function diam(rips::AbstractRipsFiltration{<:Any, T}, vertices) where T
    n = length(vertices)
    res = typemin(T)
    for i in 1:n, j in i+1:n
        d = dist(rips, vertices[j], vertices[i])
        ismissing(d) && return missing
        res = ifelse(d < res, res, d)
    end
    return ifelse(res > threshold(rips), missing, res)
end

@propagate_inbounds function diam(rips::AbstractRipsFiltration, sx::AbstractSimplex, us, v)
    res = diam(sx)
    for u in us
        # Even though this looks like a tight loop, v changes way more often than us, so
        # this is the faster order of indexing by u and v.
        d = dist(rips, v, u)
        if ismissing(d) || d > threshold(rips)
            return missing
        else
            res = ifelse(res > d, res, d)
        end
    end
    return res
end

edges(rips::AbstractRipsFiltration) = edges(dist(rips), threshold(rips), edge_type(rips))
simplex_type(rips::AbstractRipsFiltration{I, T}, dim) where {I, T} = Simplex{dim, T, I}

"""
    dist(::AbstractRipsFiltration, u, v)
    dist(::AbstractRipsFiltration)

Return the distance between vertices `u` and `v`. If the distance is higher than the
threshold, return `missing` instead. If `u` and `v` are not given, return the distance
matrix.
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
    Rips{I, T} <: AbstractRipsFiltration{I, T}

This type represents a filtration of Vietoris-Rips complexes.
Diagonal items are treated as vertex birth times.

# Constructors

* `Rips(distance_matrix; threshold=default_rips_threshold(dist))`
* `Rips(points; metric=Euclidean(), threshold=default_rips_threshold(dist))`
* `Rips{I}(args...)`: `I` sets the size of integer used to represent simplices.
"""
struct Rips{I, T, A<:AbstractMatrix{T}} <: AbstractRipsFiltration{I, T}
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
    return Rips{I}(distances(metric, points); kwargs...)
end
function Rips(dist; kwargs...)
    return Rips{Int}(dist; kwargs...)
end

@propagate_inbounds function dist(rips::Rips{<:Any, T}, i, j) where T
    return ifelse(i == j, zero(T), rips.dist[i, j])
end
dist(rips::Rips) = rips.dist

n_vertices(rips::Rips) = size(rips.dist, 1)
threshold(rips::Rips) = rips.threshold
birth(rips::Rips, i) = rips.dist[i, i]

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
struct SparseRips{I, T, A<:AbstractSparseMatrix{T}} <: AbstractRipsFiltration{I, T}
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
dist(rips::SparseRips) = rips.dist

n_vertices(rips::SparseRips) = size(rips.dist, 1)
threshold(rips::SparseRips) = rips.threshold
birth(rips::SparseRips, i) = rips.dist[i, i]
simplex_type(rips::SparseRips{I, T}, dim) where {I, T} = Simplex{dim, T, I}
