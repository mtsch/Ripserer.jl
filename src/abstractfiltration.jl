"""
    AbstractFiltration

A filtration is used to find the edges in filtration and to determine diameters of
simplices.

# Interface

* [`n_vertices(::AbstractFiltration)`](@ref)
* [`edges(::AbstractFiltration)`](@ref)
* [`diam(::AbstractFiltration, vertices)`](@ref)
* [`diam(::AbstractFiltration, ::AbstractSimplex, ::Any, ::Any)`](@ref) - only used when
  [`simplex_type`](@ref) is an [`IndexedSimplex`](@ref).
* [`simplex_type(::AbstractFiltration, dim)`](@ref)
* [`birth(::AbstractFiltration, v)`](@ref) - optional, defaults to returning `zero(T)`.
* [`threshold(::AbstractFiltration)`](@ref) - optional, defaults to returning `Inf`.
* [`postprocess_interval(::AbstractFiltration, ::Any)`](@ref) - optional
  postprocessing function that is applied to each interval in resulting persistence diagram.
"""
abstract type AbstractFiltration end

function Base.show(io::IO, flt::AbstractFiltration)
    print(io, typeof(flt), "(n_vertices=$(n_vertices(flt)))")
end

"""
    simplex_type(::AbstractFiltration, d)

Return the `d`-dimensional simplex type in the filtration.
"""
simplex_type(::AbstractFiltration, dim)

vertex_type(flt::AbstractFiltration) = simplex_type(flt, 0)
edge_type(flt::AbstractFiltration) = simplex_type(flt, 1)

"""
    n_vertices(filtration::AbstractFiltration)

Return the number of vertices in `filtration`.
"""
n_vertices(::AbstractFiltration)

"""
    edges(filtration::AbstractFiltration)

Get edges in distance matrix in `filtration`, sorted by decresing length and increasing
combinatorial index. Edges should be of type [`simplex_type`](@ref)(filtration, 1)`.
"""
edges(::AbstractFiltration)

"""
     simplex(::AbstractFiltration, ::Val{D}, vertices, sign=1)

Return `D`-simplex constructed from `vertices` with sign equal to `sign`. Return `nothing`
if simplex is not in filtration. This function is safe to call with vertices that are out of
order. Default implementation sorts `vertices` and calls [`unsafe_simplex`](@ref).
"""
function simplex(flt::AbstractFiltration, ::Val{D}, vertices, sign=1) where D
    vxs = TupleTools.sort(Tuple(vertices), rev=true)
    if allunique(vxs) && all(x -> x > 0, vxs)
        return unsafe_simplex(flt, Val(D), vxs, sign)
    else
        throw(ArgumentError("invalid vertices $(vertices)"))
    end
end

"""
    unsafe_simplex(::AbstractFiltration, ::Val{D}, vertices, sign=1)

Return `D`-simplex constructed from `vertices` with sign equal to `sign`. Return `nothing`
if simplex is not in filtration. The unsafe in the name implies that it's up to the caller
to ensure vertices are sorted.

# See also

* [`simplex`](@ref)
* [`unsafe_cofacet`](@ref)
"""
unsafe_simplex(::AbstractFiltration, ::Val, vertices, sign)

"""
    unsafe_cofacet(::AbstractFiltration, simplex, cofacet_vertices, new_vertex, sign=1)

Return cofacet of `simplex` with vertices equal to `cofacet_vertices`. `new_vertex` is the
vertex that was added to construct the cofacet. The unsafe in the name implies that it's up
to the caller to ensure vertices are sorted.

Default implementation uses [`unsafe_simplex`](@ref).
"""
function unsafe_cofacet(
    flt::AbstractFiltration, ::AbstractSimplex{D}, vertices, v, sign=1
) where D
    return unsafe_simplex(flt, Val(D + 1), vertices, sign)
end

"""
    birth(::AbstractFiltration, v)

Get the birth time of vertex `v` in filtration. Defaults to 0.
"""
birth(::AbstractFiltration, _) = false # false is a strong zero.

"""
    threshold(::AbstractFiltration)

Get the threshold of filtration. This is the maximum diameter a simplex in the filtration
can have. Used only for placing the infinity line in plotting. Defaults to `missing`.
"""
threshold(::AbstractFiltration) = Inf

"""
    postprocess_diagram(::AbstractFiltration, interval)

This function is called on each resulting persistence interval. The default implementation
does nothing.
"""
postprocess_interval(::AbstractFiltration, interval) = interval
