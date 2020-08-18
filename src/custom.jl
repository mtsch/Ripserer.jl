"""
    Custom{T} <: AbstractFiltration{Int, T}

Build a custom filtration by specifying simplices and their birth times.

The list of simplices is corrected to form a valid filtration; birth times are corrected
so a simplex is never born before its faces and missing simplices are added.

See example below for construction. Note how the unlisted 0-simplices were added with birth
times equal to the lowest between their cofaces. The order in which simplices are given does
not matter.

Note that this is not terribly efficient as Ripserer's implicit algorithm is not optimized
for these kinds of filtrations. It can still be useful for experimentation or small
filtrations.

# Example

```jldoctest
julia> flt = Custom([
    (1,) => 0,
    (4,) => 0,
    (1, 2) => 1,
    (1, 3) => 2,
    (1, 4) => 3,
    (2, 3) => 4,
    (2, 4) => 5,
    (3, 4) => 6,
    (1, 2, 3) => 7,
    (1, 2, 4) => 8,
    (1, 3, 4) => 9,
]; threshold=8)
Custom{Int64, Int64}(n_vertices=4)

julia> flt[0] # Can be indexed with dimension to list simplices
4-element Array{Simplex{0,Int64,Int64},1}:
 +Simplex{0}([4], 0)
 +Simplex{0}([2], 1)
 +Simplex{0}([3], 2)
 +Simplex{0}([1], 0)

julia> ripserer(flt)[1]
4-element 0-dimensional PersistenceDiagram:
 [0.0, 3.0)
 [1.0, ∞)

julia> ripserer(flt)[2]
3-element 1-dimensional PersistenceDiagram:
 [4.0, 7.0)
 [5.0, 8.0)
 [6.0, ∞)

"""
struct Custom{T} <: AbstractFiltration{Int, T}
    dicts::Vector{Dict{Int, T}}
    n_vertices::Int
    threshold::T
end

function Custom{T}(simplices, dim; threshold=nothing) where {T}
    dicts = [Dict{Int, T}() for _ in 0:dim]
    thresh = typemin(T)
    for (vertices, birth) in simplices
        vertices = TupleTools.sort(vertices, rev=true)
        idx = index(vertices)
        dim = length(vertices) - 1
        d_dict = dicts[dim + 1]
        d_dict[idx] = min(T(birth), get(d_dict, idx, typemax(T)))
        thresh = max(thresh, birth)

        # Propagate birth time to faces.
        for d in dim-1:-1:0
            d_dict = dicts[d + 1]
            for vs in IterTools.subsets(vertices, Val(d + 1))
                idx = index(vs)
                d_dict[idx] = min(T(birth), get(d_dict, idx, typemax(T)))
            end
        end
    end
    Custom{T}(
        dicts, maximum(keys(dicts[1])), !isnothing(threshold) ? T(threshold) : thresh
    )
end

function Custom(simplices; kwargs...)
    T = Union{}
    dim = 0
    for (vertices, birth) in simplices
        dim = max(dim, length(vertices) - 1)
        T = promote_type(T, typeof(birth))
    end
    return Custom{T}(simplices, dim; kwargs...)
end

Base.getindex(cf::Custom, d) = cf[Val(d)]
# Val for type stability with edges and column assembly.
function Base.getindex(cf::Custom, ::Val{D}) where D
    if D ≤ dim(cf)
        return [
            simplex_type(cf, D)(i, b) for (i, b) in cf.dicts[D + 1] if b < threshold(cf)
        ]
    else
        return simplex_type(cf, D)[]
    end
end

function unsafe_simplex(cf::Custom, ::Val{D}, vertices, sign) where D
    if D > dim(cf)
        return nothing
    else
        idx = index(vertices)
        birth = get(cf.dicts[D + 1], idx, nothing)
        if isnothing(birth) || birth > threshold(cf)
            return nothing
        else
            return simplex_type(cf, D)(sign * idx, birth)
        end
    end
end

dim(cf::Custom) = length(cf.dicts) - 1
simplex_type(::Type{Custom{T}}, D) where T = Simplex{D, T, Int}
n_vertices(cf::Custom) = cf.n_vertices
birth(cf::Custom, v) = cf.dicts[1][v]
birth(cf::Custom) = [cf.dicts[1][i] for i in 1:n_vertices(cf)]
edges(cf::Custom) = cf[Val(1)]
columns_to_reduce(cf::Custom, prev) = cf[Val(dim(eltype(prev)) + 1)]
threshold(cf::Custom) = cf.threshold
