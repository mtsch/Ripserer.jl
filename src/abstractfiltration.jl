"""
    AbstractFiltration{T, S<:AbstractSimplex}

A filtration is used to find the edges in filtration and to determine diameters of
simplices.

`T` is the distance type, accessible by `dist_type` and `S` is the simplex type.

# Interface

* [`n_vertices(::AbstractFiltration)`](@ref)
* [`edges(::AbstractFiltration)`](@ref)
* [`diam(::AbstractFiltration, vertices)`](@ref)
* [`diam(::AbstractFiltration, ::AbstractSimplex, ::Any, ::Any)`](@ref) - only used when
  [`simplex_type`](@ref) is an [`IndexedSimplex`](@ref).
* [`birth(::AbstractFiltration, v)`](@ref) - optional, defaults to returning `zero(T)`.
* [`threshold(::AbstractFiltration)`](@ref) - optional, defaults to returning `missing`.
* [`simplex_type(::AbstractFiltration, dim)`](@ref)
"""
abstract type AbstractFiltration{T, S<:AbstractSimplex} end

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
dist_type(::AbstractFiltration{T}) where T = T

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
    diam(::AbstractFiltration, vertices)

Get the diameter of a simplex with the vertex set `vertices`. Should return `missing` if
`vertices` do not form a valid simplex.
"""
diam(::AbstractFiltration, ::Any)

"""
    diam(::AbstractFiltration, simplex, vertices, new_vertex)

Get the diameter of coface of a `Simplex` that is formed by adding `new_vertex` to
`vertices`. Should return `missing` if new simplex is not valid.

This functions is used with the [`coboundary`](@ref) function for [`IndexedSimplex`](@ref)
"""
diam(::AbstractFiltration, ::AbstractSimplex, ::Any, ::Any)

"""
    birth(::AbstractFiltration, v)

Get the birth time of vertex `v` in filtration. Defaults to 0.
"""
birth(::AbstractFiltration{T}, _) where T = zero(T)

"""
    threshold(::AbstractFiltration)

Get the threshold of filtration. This is the maximum diameter a simplex in the filtration
can have. Used only for placing the infinity line in plotting.
"""
threshold(::AbstractFiltration) = missing
