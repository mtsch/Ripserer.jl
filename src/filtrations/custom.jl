"""
    abstract type AbstractCustomFiltration{I, T} <: AbstractFiltration{I, T}

This abstract type is for filtrations that have all simplices stored in `Dict`s. The dicts
should be accessible by the function [`simplex_dicts`](@ref) and should be a vector of
`Dict{I, T}`. A custom filtration should also have [`adjacency_matrix`](@ref) defined. This
matrix is only used as an adjacency matrix. Its values are ignored.
"""
abstract type AbstractCustomFiltration{I,T} <: AbstractFiltration{I,T} end

"""
    simplex_dicts(::AbstractCustomFiltration)

Get the dictionaries used to get simplex birth times. Should return a `Vector` of `Dict{I,
T}` that maps a simplex index to its birth time. The first element of this `Vector`
corresponds to vertices, second to 1-simplices etc.
"""
simplex_dicts

Base.getindex(cf::AbstractCustomFiltration, d) = cf[Val(d)]
# Val for type stability with edges and column assembly.
function Base.getindex(cf::AbstractCustomFiltration, ::Val{D}) where {D}
    if D ≤ dim(cf)
        return [
            simplex_type(cf, D)(i, b) for
            (i, b) in simplex_dicts(cf)[D + 1] if b ≤ threshold(cf)
        ]
    else
        return simplex_type(cf, D)[]
    end
end

function unsafe_simplex(
    ::Type{S}, cf::AbstractCustomFiltration, vertices
) where {D,S<:Simplex{D}}
    if D > dim(cf)
        return nothing
    else
        idx = index(vertices)
        birth = get(simplex_dicts(cf)[D + 1], idx, nothing)
        if isnothing(birth) || birth > threshold(cf)
            return nothing
        else
            return simplex_type(cf, D)(idx, birth)
        end
    end
end

dim(cf::AbstractCustomFiltration) = length(simplex_dicts(cf)) - 1
simplex_type(::Type{<:AbstractCustomFiltration{I,T}}, D) where {I,T} = Simplex{D,T,I}
births(cf::AbstractCustomFiltration) = [simplex_dicts(cf)[1][i] for i in 1:nv(cf)]
edges(cf::AbstractCustomFiltration) = cf[Val(1)]
columns_to_reduce(cf::AbstractCustomFiltration, prev) = cf[Val(dim(eltype(prev)) + 1)]
nv(cf::AbstractCustomFiltration) = size(adjacency_matrix(cf), 1)
coboundary(cf::AbstractCustomFiltration, sx::Simplex) = SparseCoboundary{true}(cf, sx)

"""
    Custom{I, T} <: AbstractCustomFiltration{I, T}

Build a custom filtration by specifying simplices and their birth times.

The list of simplices is corrected to form a valid filtration; birth times are corrected
so a simplex is never born before its faces and missing simplices are added.

See the examples below for construction. Note how the unlisted 0-simplices were added with
birth times equal to the lowest between their cofaces. The order in which simplices are
given does not matter.

To create your own types of custom filtrations, subtype [`AbstractCustomFiltration`](@ref).

# Examples

```jldoctest
julia> flt = Custom([(1,) => 0, (4,) => 0, (1, 2) => 1, (1, 3) => 2, (1, 4) => 3, (2, 3) => 4, (2, 4) => 5, (3, 4) => 6, (1, 2, 3) => 7, (1, 2, 4) => 8, (1, 3, 4) => 9]; threshold=8)
Custom{Int64, Int64}(nv=4)

julia> flt[0] # Can be indexed with dimension to list simplices
4-element Vector{Simplex{0, Int64, Int64}}:
 +Simplex{0}((4,), 0)
 +Simplex{0}((2,), 1)
 +Simplex{0}((3,), 2)
 +Simplex{0}((1,), 0)

julia> ripserer(flt)[1]
2-element 0-dimensional PersistenceDiagram:
 [0.0, 3.0)
 [0.0, ∞)

julia> ripserer(flt)[2]
3-element 1-dimensional PersistenceDiagram:
 [5.0, 8.0)
 [4.0, 7.0)
 [6.0, ∞)

```
"""
struct Custom{I,T} <: AbstractCustomFiltration{I,T}
    adj::SparseMatrixCSC{Bool,Int} # adjacency matrix for sparse coboundary
    dicts::Vector{Dict{I,T}}
    threshold::T
end

function insert_simplex!(
    dicts::Vector{Dict{I,T}}, vertices::NTuple{N,I}, birth, threshold
) where {N,T,I}
    if birth > threshold
        return nothing
    else
        idx = index(vertices)
        dim = N - 1
        d_dict = dicts[dim + 1]
        d_dict[idx] = min(birth, get(d_dict, idx, typemax(T)))
        _insert_simplex_facets!(dicts, vertices, birth, Val(N - 1))
    end
end
@inline _insert_simplex_facets!(_, _, _, ::Val{0}) = nothing
@inline function _insert_simplex_facets!(dicts, vertices, birth, ::Val{N}) where {N}
    dim = N - 1
    d_dict = dicts[dim + 1]
    T = valtype(d_dict)
    for vs in IterTools.subsets(vertices, Val(N))
        idx = index(vs)
        d_dict[idx] = min(birth, get(d_dict, idx, typemax(T)))
    end
    return _insert_simplex_facets!(dicts, vertices, birth, Val(N - 1))
end

function _adjacency_matrix(dicts)
    nv = maximum(keys(dicts[1]))
    adj_is = Int[]
    adj_js = Int[]
    adj_vs = Bool[]
    for (index, _) in dicts[2]
        u, v = _vertices(index, Val(2))
        append!(adj_is, (u, v))
        append!(adj_js, (v, u))
        append!(adj_vs, (true, true))
    end
    return sparse(adj_is, adj_js, adj_vs, nv, nv)
end

function Custom{I,T}(simplices, dim_max::Int, threshold::T) where {I,T}
    dicts = [Dict{I,T}() for _ in 0:dim_max]
    threshold = T(threshold)

    for (_vertices, birth) in simplices
        vertices = I.(TupleTools.sort(_vertices; rev=true))
        if birth ≤ threshold
            insert_simplex!(dicts, vertices, birth, threshold)
        end
    end
    adj = _adjacency_matrix(dicts)

    return Custom{I,T}(adj, dicts, threshold)
end

# TODO: hot mess. Simplex sorting and index conversion could be done here.
function Custom{I}(simplices; threshold=nothing, verbose=false) where {I}
    # Promote birth types and find dim and threshold.
    T = Union{}
    dim = 0
    thresh = typemin(simplices[1][2])
    largest_simplex = ()
    for (vertices, birth) in simplices
        dim = max(dim, length(vertices) - 1)
        T = promote_type(T, typeof(birth))
        if isnothing(threshold)
            thresh = max(thresh, birth)
        end
        if length(vertices) ≥ length(largest_simplex)
            largest_simplex = max(TupleTools.sort(vertices; rev=true), largest_simplex)
        end
    end
    if isnothing(threshold)
        threshold = thresh
    end
    index_overflow_check(I.(largest_simplex))
    return Custom{I,T}(simplices, dim, T(threshold))
end

Custom(args...; kwargs...) = Custom{Int}(args...; kwargs...)

simplex_dicts(cf::Custom) = cf.dicts
adjacency_matrix(cf::Custom) = cf.adj
threshold(cf::Custom) = cf.threshold
