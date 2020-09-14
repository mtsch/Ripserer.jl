"""
    AbstractRipsFiltration{I<:Signed, T} <: AbstractFiltration{I, T}

An abstract Vietoris-Rips filtration. Its subtypes can only overload
[`adjacency_matrix`](@ref) and get default implementations for the rest of the filtration
interface.

# Example

```jldoctest
julia> struct MyRips <: Ripserer.AbstractRipsFiltration{Int, Float16} end

julia> Ripserer.adjacency_matrix(::MyRips) = [0 1 1; 1 0 1; 1 1 0]

julia> ripserer(MyRips())
2-element Array{PersistenceDiagrams.PersistenceDiagram,1}:
 3-element 0-dimensional PersistenceDiagram
 0-element 1-dimensional PersistenceDiagram

```
"""
abstract type AbstractRipsFiltration{I<:Signed, T} <: AbstractFiltration{I, T} end

nv(rips::AbstractRipsFiltration) = size(adjacency_matrix(rips), 1)

function births(rips::AbstractRipsFiltration)
    adj = adjacency_matrix(rips)
    return view(adj, diagind(adj))
end

simplex_type(::Type{<:AbstractRipsFiltration{I, T}}, D) where {I, T} = Simplex{D, T, I}

function edges(rips::AbstractRipsFiltration)
    if issparse(adjacency_matrix(rips))
        _sparse_edges(rips)
    else
        _full_edges(rips)
    end
end

function _full_edges(rips::AbstractRipsFiltration)
    result = edge_type(rips)[]
    @inbounds for u in 1:nv(rips), v in u+1:nv(rips)
        sx = unsafe_simplex(rips, Val(1), (v, u), 1)
        !isnothing(sx) && push!(result, sx)
    end
    return result
end

function _sparse_edges(rips::AbstractRipsFiltration)
    adj = adjacency_matrix(rips)
    result = edge_type(rips)[]
    rows = rowvals(adj)
    for u in 1:nv(rips)
        for i in nzrange(adj, u)
            v = rows[i]
            if v > u
                sx = unsafe_simplex(rips, Val(1), (v, u), 1)
                !isnothing(sx) && push!(result, sx)
            end
        end
    end
    return result
end

@inline @propagate_inbounds function unsafe_simplex(
    ::Type{S}, rips::AbstractRipsFiltration{I, T}, vertices, sign
) where {I, T, S<:Simplex}
    if dim(S) == 0
        return S(I(sign) * vertices[1], births(rips)[vertices[1]])
    else
        adj = adjacency_matrix(rips)
        n = length(vertices)
        diameter = typemin(T)
        for i in 1:n, j in i+1:n
            d = adj[vertices[j], vertices[i]]
            (iszero(d) || d > threshold(rips)) && return nothing
            diameter = ifelse(d > diameter, d, diameter)
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
    adj = adjacency_matrix(rips)
    for v in cofacet_vertices
        v == new_vertex && continue
        # Even though this looks like a tight loop, v changes way more often than us, so
        # this is the faster order of indexing by new_vertex and v.
        d = adj[new_vertex, v]
        (iszero(d) || d > threshold(rips)) && return nothing
        diameter = ifelse(d > diameter, d, diameter)
    end
    return S(I(sign) * index(cofacet_vertices), diameter)
end

# This definition is used in SparseCoboundary
@inline @propagate_inbounds function unsafe_cofacet(
    ::Type{S},
    rips::AbstractRipsFiltration{I},
    simplex,
    cofacet_vertices,
    _,
    sign,
    new_edges,
) where {I, S<:Simplex}
    diameter = birth(simplex)
    for e in new_edges
        e > threshold(rips) && return nothing
        diameter = ifelse(e > diameter, e, diameter)
    end
    return S(I(sign) * index(cofacet_vertices), diameter)
end

function coboundary(
    rips::AbstractRipsFiltration, sx::AbstractSimplex, ::Val{A}=Val(true)
) where A
    if adjacency_matrix(rips) isa SparseMatrixCSC
        return SparseCoboundary{A}(rips, sx)
    else
        return Coboundary{A}(rips, sx)
    end
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
    radius(points[, metric=Euclidean(1e-12)])

Calculate the radius of the space. This is used for default `thresholds`.
"""
function radius(dists::AbstractMatrix)
    return minimum(maximum(abs, dists[:, i]) for i in 1:size(dists, 1))
end
function radius(dists::SparseMatrixCSC)
    return maximum(dists)
end
function radius(points, metric=Euclidean(1e-12))
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

function _check_distance_matrix(dist::SparseMatrixCSC)
    n = LinearAlgebra.checksquare(dist)
    vals = nonzeros(dist)
    rows = rowvals(dist)
    for i in 1:n
        vertex_birth = dist[i, i]
        for j in nzrange(dist, i)
            i == rows[j] && continue
            edge_birth = vals[j]
            iszero(edge_birth) && throw(ArgumentError(
                "zero edges in input matrix"
            ))
            edge_birth < vertex_birth && throw(ArgumentError(
                "edges with birth value lower than vertex births in input matrix"
            ))
            edge_birth ≠ dist[i, rows[j]] && throw(ArgumentError(
                "input matrix not symmetric"
            ))
        end
    end
end

function _check_distance_matrix(dist)
    n = LinearAlgebra.checksquare(dist)
    for i in 1:n
        for j in i+1:n
            vertex_birth = max(dist[j, j], dist[i, i])
            edge_birth = dist[j, i]
            iszero(edge_birth) && throw(ArgumentError(
                "zero edges in input matrix"
            ))
            edge_birth < vertex_birth && throw(ArgumentError(
                "edges with birth value lower than vertex births in input matrix"
            ))
            edge_birth ≠ dist[i, j] && throw(ArgumentError(
                "input matrix not symmetric"
            ))
        end
    end
end

"""
    Rips{I, T} <: AbstractRipsFiltration{I, T}

This type represents a filtration of Vietoris-Rips complexes.

Diagonal items in the input matrix are treated as vertex birth times.

Zero values are not allowed due to how sparse matrices work in Julia. If you need zero birth
times, try offseting all values by a constant.

Threshold defaults to the radius of the input space. When using low `threshold`s, consider
using the `sparse=true` keyword argument. It will give the same result, but may be much
faster.

# Constructors

* `Rips(distance_matrix; threshold=nothing)`
* `Rips(points; metric=Euclidean(1e-12), threshold=nothing)`
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
    adj::A
    threshold::T
end

function Rips{I}(
    dists::AbstractMatrix{T}; threshold=nothing, sparse=false,
) where {I, T}
    _check_distance_matrix(dists)
    thresh = isnothing(threshold) ? radius(dists) : T(threshold)
    if sparse
        adj = SparseArrays.sparse(dists)
    else
        adj = copy(dists)
    end
    if issparse(adj)
        adj isa SparseMatrixCSC || throw(ArgumentError(
            "only SparseMatrixCSC sparse matrices supported"
        ))
        SparseArrays.fkeep!(adj, (_, _, v) -> v ≤ thresh)
    end
    return Rips{I, T, typeof(adj)}(adj, thresh)
end
function Rips{I}(points::AbstractVector; metric=Euclidean(1e-12), kwargs...) where I
    return Rips{I}(distances(metric, unique(points)); kwargs...)
end
function Rips(dist_or_points; kwargs...)
    return Rips{Int}(dist_or_points; kwargs...)
end

function Base.show(io::IO, rips::Rips{I, T}) where {I, T}
    print(io, "Rips{$I, $T}(nv=$(nv(rips)), sparse=$(issparse(rips.adj)))")
end

threshold(rips::Rips) = rips.threshold
adjacency_matrix(rips::Rips) = rips.adj

@deprecate SparseRips(args...; kwargs...) Rips(args...; sparse=true, kwargs...)
