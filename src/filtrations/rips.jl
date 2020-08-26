"""
    AbstractRipsFiltration{I<:Signed, T} <: AbstractFiltration{I, T}

An abstract Vietoris-Rips filtration. Its subtypes can only overload [`dist`](@ref) and get
default implementations for the rest of the filtration interface.

# Example

```jldoctest
julia> struct MyRips <: Ripserer.AbstractRipsFiltration{Int, Float16} end

julia> Ripserer.dist(::MyRips) = [0 1 1; 1 0 1; 1 1 0]

julia> ripserer(MyRips())
2-element Array{PersistenceDiagrams.PersistenceDiagram,1}:
 3-element 0-dimensional PersistenceDiagram
 0-element 1-dimensional PersistenceDiagram

```
"""
abstract type AbstractRipsFiltration{I<:Signed, T} <: AbstractFiltration{I, T} end

"""
    dist(::AbstractRipsFiltration, u, v)
    dist(::AbstractRipsFiltration)

Return the distance between vertices `u` and `v`. If the distance is somehow invalid, it may
return `missing` instead. If `u` and `v` are not given, return the distance matrix.

# Examples

```jldoctest
julia> flt = Rips([1 1 2; 1 0 1; 2 1 0]);

julia> Ripserer.dist(flt, 1, 3)
2

julia> Ripserer.dist(flt)
3×3 Array{Int64,2}:
 1  1  2
 1  0  1
 2  1  0

```
"""
dist(rips::AbstractRipsFiltration, i, j) = dist(rips)[i, j]

nv(rips::AbstractRipsFiltration) = size(dist(rips), 1)
birth(rips::AbstractRipsFiltration) = diag(dist(rips))
simplex_type(::Type{<:AbstractRipsFiltration{I, T}}, D) where {I, T} = Simplex{D, T, I}
adjacency_matrix(rips::AbstractRipsFiltration) = dist(rips)

function edges(rips::AbstractRipsFiltration)
    if issparse(dist(rips))
        _sparse_edges(rips)
    else
        _full_edges(rips)
    end
end

function _full_edges(rips::AbstractRipsFiltration)
    result = edge_type(rips)[]
    @inbounds for u in 1:size(dist(rips), 1), v in u+1:size(dist(rips), 1)
        sx = unsafe_simplex(rips, Val(1), (v, u), 1)
        !isnothing(sx) && push!(result, sx)
    end
    return result
end

function _sparse_edges(rips::AbstractRipsFiltration)
    result = edge_type(rips)[]
    rows = rowvals(rips.dist)
    vals = nonzeros(rips.dist)
    for u in 1:size(rips.dist, 1)
        for i in nzrange(dist(rips), u)
            v = rows[i]
            if v > u
                sx = unsafe_simplex(rips, Val(1), (v, u), 1)
                !isnothing(sx) && push!(result, sx)
            end
        end
    end
    return result
end

"""
    to_matrix(points)

Convert collection of d-dimensional points to d×n matrix
"""
function to_matrix(points)
    dim = length(points[1])
    T = eltype(points[1])
    n = length(points)
    result = zeros(T, (dim, n))

    for (i, p) in enumerate(points)
        length(p) == dim || throw(ArgumentError("points must have the same length"))
        result[:, i] .= p
    end

    return result
end

"""
    distances(metric, points)

Return distance matrix calculated from `points` with `metric`.
"""
function distances(metric, points)
    points_mat = to_matrix(points)
    dists = pairwise(metric, points_mat, dims=2)
    return dists
end

"""
    radius(dists)
    radius(points[, metric=Euclidean()])

Calculate the radius of the space. This is used for default `thresholds`.
"""
function radius(dists::AbstractMatrix{T}) where T
    return minimum(maximum(abs, dists[:, i]) for i in 1:size(dists, 1))
end
function radius(points, metric=Euclidean())
    radius = Inf
    for p in points
        p_max = 0.0
        for q in points
            p == q && continue
            p_max = max(p_max, metric(SVector(p), SVector(q)))
        end
        radius = min(p_max, radius)
    end
    return radius
end

@inline @propagate_inbounds function unsafe_simplex(
    ::Type{S}, rips::AbstractRipsFiltration{I, T}, vertices, sign
) where {I, T, S<:Simplex}
    if dim(S) == 0
        return S(I(sign) * vertices[1], birth(rips, vertices[1]))
    else
        n = length(vertices)
        diameter = typemin(T)
        for i in 1:n, j in i+1:n
            d = dist(rips, vertices[j], vertices[i])
            if ismissing(d) || d > threshold(rips)
                return nothing
            else
                _d::T = d
                diameter = ifelse(_d > diameter, _d, diameter)
            end
        end
        return S(I(sign) * index(vertices), diameter)
    end
end

@inline @propagate_inbounds function unsafe_cofacet(
    ::Type{S},
    rips::AbstractRipsFiltration{I, T},
    simplex,
    cofacet_vertices,
    new_vertex,
    sign,
) where {I, T, S<:Simplex}
    diameter = birth(simplex)
    for v in cofacet_vertices
        v == new_vertex && continue
        # Even though this looks like a tight loop, v changes way more often than us, so
        # this is the faster order of indexing by new_vertex and v.
        d = dist(rips, new_vertex, v)
        if ismissing(d) || d > threshold(rips)
            return nothing
        else
            _d::T = d
            diameter = ifelse(_d > diameter, _d, diameter)
        end
    end
    return S(I(sign) * index(cofacet_vertices), diameter)
end

# This definition is used in SparseCoboundary
@inline @propagate_inbounds function unsafe_cofacet(
    ::Type{S},
    rips::AbstractRipsFiltration{I},
    sx,
    cofacet_vertices,
    _,
    sign,
    new_edges,
) where {I, S<:Simplex}
    new_diam = birth(sx)
    for e in new_edges
        e > threshold(rips) && return nothing
        new_diam = ifelse(e > new_diam, e, new_diam)
    end
    new_index = index(cofacet_vertices)
    return S(I(sign) * new_index, new_diam)
end

function coboundary(
    rips::AbstractRipsFiltration, sx::AbstractSimplex, ::Val{A}=Val(true)
) where A
    if dist(rips) isa SparseMatrixCSC
        return SparseCoboundary{A}(rips, sx)
    else
        return Coboundary{A}(rips, sx)
    end
end

"""
    Rips{I, T} <: AbstractRipsFiltration{I, T}

This type represents a filtration of Vietoris-Rips complexes.
Diagonal items are treated as vertex birth times.

Threshold defaults to the radius of the input space. When using low `threshold`s, consider
using [`SparseRips`](@ref) instead. It will give the same result, but may be much faster.

# Constructors

* `Rips(distance_matrix; threshold=nothing)`
* `Rips(points; metric=Euclidean(), threshold=nothing)`
* `Rips{I}(args...)`: `I` sets the size of integer used to represent simplices.

# Examples

```jldoctest
julia> data = [(sin(t), cos(t)) for t in range(0, 2π, length=101)][1:end-1];

julia> ripserer(Rips(data))
2-element Array{PersistenceDiagrams.PersistenceDiagram,1}:
 100-element 0-dimensional PersistenceDiagram
 1-element 1-dimensional PersistenceDiagram

julia> ripserer(Rips(data, threshold=1.7))[2]
1-element 1-dimensional PersistenceDiagram:
 [0.0628, ∞)

julia> using Distances

julia> ripserer(Rips(data, metric=Cityblock()))[2]
1-element 1-dimensional PersistenceDiagram:
 [0.0888, 2.0)

```
"""
struct Rips{I, T, A<:AbstractMatrix{T}} <: AbstractRipsFiltration{I, T}
    dist::A
    threshold::T
end

function Rips{I}(
    dist::AbstractMatrix{T}; threshold=nothing
) where {I, T}
    issymmetric(dist) ||
        throw(ArgumentError("`dist` must be symmetric"))
    !issparse(dist) ||
        throw(ArgumentError("`dist` is sparse. Use `SparseRips` instead"))
    thresh = isnothing(threshold) ? radius(dist) : T(threshold)
    return Rips{I, T, typeof(dist)}(dist, thresh)
end
function Rips{I}(points::AbstractVector; metric=Euclidean(), kwargs...) where I
    return Rips{I}(distances(metric, points); kwargs...)
end
function Rips(dist; kwargs...)
    return Rips{Int}(dist; kwargs...)
end

@propagate_inbounds function dist(rips::Rips{<:Any, T}, i, j) where T
    return rips.dist[i, j]
end
dist(rips::Rips) = rips.dist

threshold(rips::Rips) = rips.threshold
birth(rips::Rips, i) = rips.dist[i, i]

# sparse rips filtration ================================================================= #
"""
    SparseRips{I, T} <: AbstractRipsFiltration{T, Simplex}

This type represents a sparse filtration of Vietoris-Rips complexes.
The distance matrix will be converted to a sparse matrix with all values greater than
threshold deleted. Off-diagonal zeros in the matrix are treated as `missing`. Diagonal items
are treated as vertex birth times.

Use this filtration type with low `threshold`s.

# Constructor

* `SparseRips{I}(distance_matrix; threshold=nothing)`
* `SparseRips(distance_matrix; threshold=nothing)`: `I` sets the integer size used to
  represent simplices.

# Examples

```jldoctest
julia> data = [(sin(t), cos(t)) for t in range(0, 2π, length=101)][1:end-1];

julia> ripserer(SparseRips(data))
2-element Array{PersistenceDiagrams.PersistenceDiagram,1}:
 100-element 0-dimensional PersistenceDiagram
 1-element 1-dimensional PersistenceDiagram

julia> ripserer(SparseRips(data, threshold=1.7))[2]
1-element 1-dimensional PersistenceDiagram:
 [0.0628, ∞)

julia> using Distances

julia> ripserer(SparseRips(data, metric=Cityblock()))[2]
1-element 1-dimensional PersistenceDiagram:
 [0.0888, 2.0)

```
"""
struct SparseRips{I, T} <: AbstractRipsFiltration{I, T}
    dist::SparseMatrixCSC{T, Int}
    threshold::T
end

function SparseRips{I}(
    dist::AbstractMatrix{T}; threshold=nothing
) where {I, T}
    issymmetric(dist) || throw(ArgumentError("`dist` must be symmetric"))
    if isnothing(threshold)
        threshold = issparse(dist) ? maximum(dist) : radius(dist)
    end
    dists = SparseArrays.fkeep!(SparseMatrixCSC(dist), (_, _, v) -> v ≤ threshold)

    return SparseRips{I, T}(dists, threshold)
end
function SparseRips(dist; kwargs...)
    return SparseRips{Int}(dist; kwargs...)
end
function SparseRips{I}(points::AbstractVector; metric=Euclidean(), kwargs...) where I
    return SparseRips{I}(distances(metric, points); kwargs...)
end

@propagate_inbounds function dist(rips::SparseRips{<:Any, T}, i, j) where T
    res = rips.dist[i, j]
    return ifelse(i == j, res, ifelse(iszero(res), missing, res))
end
dist(rips::SparseRips) = rips.dist

threshold(rips::SparseRips) = rips.threshold
birth(rips::SparseRips, i) = rips.dist[i, i]
simplex_type(::SparseRips{I, T}, dim) where {I, T} = Simplex{dim, T, I}
