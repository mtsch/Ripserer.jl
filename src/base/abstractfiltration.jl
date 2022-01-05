"""
    AbstractFiltration{I,T}

A filtration is used to find the edges in filtration and to create simplices. An
`AbstractFiltration{I,T}`'s simplex type is expected to return simplices of type
`<:AbstractCell{_,T,I}`.

An `AbstractFiltration` constructor must accept the `verbose` keyword argument.

# Interface

* [`nv(::AbstractFiltration)`](@ref)
* [`births(::AbstractFiltration)`](@ref)
* [`vertices(::AbstractFiltration)`](@ref)
* [`edges(::AbstractFiltration)`](@ref)
* [`simplex_type(::Type{AbstractFiltration}, dim)`](@ref)
* [`simplex(::AbstractFiltration, ::Val{dim}, vertices)`](@ref)
* [`unsafe_simplex(::AbstractFiltration, ::Val{dim}, vertices)`](@ref)
* [`unsafe_cofacet`](@ref)`(::AbstractFiltration, simplex, vertices, vertex[, edges])`
* [`threshold(::AbstractFiltration)`](@ref)
* [`columns_to_reduce(::AbstractFiltration, ::Any)`](@ref)
* [`emergent_pairs(::AbstractFiltration)`](@ref)
* [`postprocess_diagram(::AbstractFiltration, ::Any)`](@ref)
"""
abstract type AbstractFiltration{I,T} end

function Base.show(io::IO, flt::AbstractFiltration{I,T}) where {I,T}
    return print(io, nameof(typeof(flt)), "{$I, $T}(nv=$(nv(flt)))")
end

"""
    simplex_type(::Type{<:AbstractFiltration}, D)
    simplex_type(::AbstractFiltration, D)

Return the `D`-dimensional simplex type in the filtration. Only the method for the type
needs to be overloaded.

# Examples

```jldoctest
julia> Ripserer.simplex_type(Rips{Int,Float64}, 1)
Simplex{1, Float64, Int64}

julia> Ripserer.simplex_type(Cubical{2,Float16}, 2)
Cube{2, Float16, 2}

```
"""
simplex_type(flt::AbstractFiltration, dim) = simplex_type(typeof(flt), dim)

vertex_type(flt::AbstractFiltration) = simplex_type(typeof(flt), 0)
edge_type(flt::AbstractFiltration) = simplex_type(typeof(flt), 1)

"""
    nv(::AbstractFiltration)

Return the number of vertices in `filtration`.

# Example

```jldoctest
julia> Ripserer.nv(Rips([1 1; 1 1]))
2

```
"""
nv(::AbstractFiltration)

"""
    edges(::AbstractFiltration)

Get edges (1-simplices) in `filtration`. Edges should be of type
[`simplex_type`](@ref)`(filtration, 1)`.

# Example

```
julia> Ripserer.edges(Rips([0 2 1; 2 0 1; 1 1 0], threshold=2))
3-element Array{Simplex{1,Int64,Int64},1}:
 +Simplex{1}([2, 1], 2)
 +Simplex{1}([3, 1], 1)
 +Simplex{1}([3, 2], 1)

```
"""
edges(::AbstractFiltration)

"""
     simplex(::AbstractFiltration, ::Val{D}, vertices)

Return `D`-simplex constructed from `vertices`. Return `nothing`
if simplex is not in filtration. This function is safe to call with vertices that are out of
order. Default implementation sorts `vertices` and calls [`unsafe_simplex`](@ref).

# Example

```
julia> simplex(Rips([0 2 1; 2 0 1; 1 1 0], threshold=2), Val(1), (1, 2))
1-dimensional Simplex(index=1, birth=2):
  +[2, 1]

```
"""
function simplex(flt::AbstractFiltration, ::Val{D}, vertices) where {D}
    if D == 0 && length(vertices) == 1 && vertices[1] > 0
        # No need to run allunique in this case.
        return unsafe_simplex(flt, Val(0), vertices)
    else
        vs = TupleTools.sort(Tuple(vertices); rev=true)
        if allunique(vs) &&
            all(x -> x > 0, vs) &&
            length(vs) == length(simplex_type(flt, D))
            return unsafe_simplex(flt, Val(D), vs)
        end
    end
    return throw(ArgumentError("invalid vertices $(vertices)"))
end

"""
    unsafe_simplex(::AbstractFiltration, ::Val{D}, vertices)

Return `D`-simplex constructed from `vertices`. Return `nothing`
if simplex is not in filtration. The unsafe in the name implies that it's up to the caller
to ensure vertices are sorted and unique.
"""
function unsafe_simplex(flt, ::Val{D}, vertices) where {D}
    return unsafe_simplex(simplex_type(typeof(flt), D), flt, vertices)
end

"""
    unsafe_cofacet(filtration, simplex, cofacet_vertices, v[, edges])
    unsafe_cofacet(::Type{S}, filtration, simplex, cofacet_vertices, v[, edges])

Return cofacet of `simplex` with vertices equal to `cofacet_vertices`. `v` is the
vertex that was added to construct the cofacet. In the case of sparse rips filtrations, an
additional argument `edges` is used. `edges` is a vector that contains the weights on edges
connecting the new vertex to old vertices. `S` is the simplex type which can be used for
dispatch.

The unsafe in the name implies that it's up to the caller to ensure vertices are sorted and
unique.

Default implementation uses [`unsafe_simplex`](@ref).
"""
@inline @propagate_inbounds function unsafe_cofacet(flt, sx, args...)
    return unsafe_cofacet(simplex_type(typeof(flt), dim(sx) + 1), flt, sx, args...)
end
function unsafe_cofacet(::Type{S}, flt, _, vertices, _, args...) where {S}
    return unsafe_simplex(flt, Val(dim(S)), vertices)::Union{S,Nothing}
end

"""
    vertices(::AbstractFiltration)

Return the vertices in filtration. Defaults to `1:n`. The `shape` and `eltype` of the result
can be anything as long as `result[result[i]] == result[i]` holds.

# Example

```jldoctest
julia> vertices(Rips([0 1 1; 1 0 1; 1 1 0]))
Base.OneTo(3)

julia> vertices(Cubical([0 1 1; 1 0 1; 1 1 0]))
3×3 CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}}:
 CartesianIndex(1, 1)  CartesianIndex(1, 2)  CartesianIndex(1, 3)
 CartesianIndex(2, 1)  CartesianIndex(2, 2)  CartesianIndex(2, 3)
 CartesianIndex(3, 1)  CartesianIndex(3, 2)  CartesianIndex(3, 3)

```
"""
vertices(flt::AbstractFiltration) = Base.OneTo(nv(flt))

"""
    births(::AbstractFiltration)

Get the birth times of vertices in filtration. Defaults to all births being 0. Must return
array of the same shape as the filtration's `vertices`.

# Examples

```jldoctest
julia> flt = Rips([1 1 2; 1 0 1; 2 1 0]);

julia> Ripserer.births(flt)
3-element view(::Vector{Int64}, 1:4:9) with eltype Int64:
 1
 0
 0

```
"""
births(flt::AbstractFiltration{<:Any,T}) where {T} = zeros(T, size(vertices(flt)))

"""
    threshold(::AbstractFiltration)

Get the threshold of filtration. This is the maximum diameter a simplex in the filtration
can have. Defaults to `Inf`.

# Examples

```jldoctest
julia> threshold(Rips([0 2 1; 2 0 1; 1 1 0]))
1

julia> threshold(Cubical([1 1 2; 3 2 1; 0 0 0]))
3

```
"""
threshold(::AbstractFiltration) = Inf

"""
    adjacency_matrix(::AbstractFiltration)

Return the adjacency matrix. For sparse filtrations, this should return a `SparseMatrixCSC`.

# Examples

```jldoctest
julia> Ripserer.adjacency_matrix(Rips([0 2 1; 2 0 1; 1 1 0]))
3×3 Matrix{Int64}:
 0  2  1
 2  0  1
 1  1  0

julia> Ripserer.adjacency_matrix(Rips([0 10 2; 10 0 1; 2 1 0]; sparse=true))
3×3 SparseArrays.SparseMatrixCSC{Int64, Int64} with 4 stored entries:
 ⋅  ⋅  2
 ⋅  ⋅  1
 2  1  ⋅

julia> Ripserer.adjacency_matrix(Custom([(2, 1) => 1, (5, 1) => 2, (3, 4) => 3]))
5×5 SparseArrays.SparseMatrixCSC{Bool, Int64} with 6 stored entries:
 ⋅  1  ⋅  ⋅  1
 1  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  1  ⋅
 ⋅  ⋅  1  ⋅  ⋅
 1  ⋅  ⋅  ⋅  ⋅

```
"""
adjacency_matrix(::AbstractFiltration)

"""
    columns_to_reduce(::AbstractFilration, prev_column_itr)

List all columns to reduce in next dimension, possibly computing it from previous
columns. Default implementation uses [`coboundary`](@ref) with the all cofacets parameter
set to `Val(false)`.

# Example

```jldoctest
julia> flt = Rips([0 1 1; 1 0 1; 1 1 0]);

julia> Ripserer.columns_to_reduce(flt, Ripserer.edges(flt)) |> collect
1-element Vector{Simplex{2, Int64, Int64}}:
 +Simplex{2}((3, 2, 1), 1)

```
"""
function columns_to_reduce(flt::AbstractFiltration, prev_column_itr)
    return Iterators.flatten(imap(σ -> coboundary(flt, σ, Val(false)), prev_column_itr))
end

"""
    emergent_pairs(::AbstractFiltration)

Return `true` if the emergent pairs optimization is to be performed. Default to returning
`true`. Should be set to `false` for a filtration type that is unable to produce
(co)boundary simplices in the correct order.
"""
emergent_pairs(::AbstractFiltration) = true

"""
    postprocess_diagram(::AbstractFiltration, diagram)

This function is called on each resulting persistence diagram after all intervals have been
computed. Defaults to not doing anything.
"""
postprocess_diagram(::AbstractFiltration, diagram) = diagram

"""
    find_apparent_pairs(::AbstractFiltration, columns, is_cohomology, verbose)

Find apparent pairs in filtration. Columns are the initial columns that are being reduced.
Should return a tuple of:

* A collection of columns that are left to reduce.
* A collection of apparent pairs `(σ, τ)`. These will be insered into the reduced matrix so
  that `matrix[τ] = [σ]`.

!!! note
    Mutating `columns` and returning it from this function is safe.

!!! note
    This can only be applied to cohomology.

!!! warning
    This is still experimental.
"""
find_apparent_pairs(::AbstractFiltration, columns, _) = columns, ()

function index_overflow_check(
    flt::AbstractFiltration,
    ::Type{F},
    dim_max,
    message="$flt overflows in $(dim_max)D. " *
            "Try using a bigger index type or lower `dim_max`",
) where {F}
    S = simplex_type(flt, dim_max + 1)
    length(S) > nv(flt) && throw(ArgumentError("$S has more than $(nv) vertices."))
    return index_overflow_check(S, F, nv(flt), message)
end

"""
    distance_matrix(::AbstractFiltration)

Return a matrix with distances between vertices of the filtration. These distances are used
to determine edge length when finding shortest represenatative cycle in
[`reconstruct_cycle`](@ref).

Defaults to all distances being 1.

# Examples

```jldoctest
julia> flt = Rips([1 1 2; 1 0 1; 2 1 0]);

julia> Ripserer.distance_matrix(flt)
3×3 Matrix{Int64}:
 1  1  2
 1  0  1
 2  1  0

julia> flt = Custom([(1,2,3,4) => 10.0, (1,2) => 3.0, (1,3) => 4.0, (2,3) => 5.0]);

julia> Ripserer.distance_matrix(flt)
4×4 Ripserer.DefaultDist{Float64}:
 0.0  1.0  1.0  1.0
 1.0  0.0  1.0  1.0
 1.0  1.0  0.0  1.0
 1.0  1.0  1.0  0.0

```
"""
distance_matrix(flt::AbstractFiltration{<:Any,T}) where {T} = DefaultDist{T}(nv(flt))

struct DefaultDist{T} <: AbstractMatrix{T}
    nv::Int
end

Base.size(dd::DefaultDist) = (dd.nv, dd.nv)
Base.getindex(dd::DefaultDist, i, j) = ifelse(i == j, zero(eltype(dd)), one(eltype(dd)))
