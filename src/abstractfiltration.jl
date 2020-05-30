"""
    AbstractFiltration{T, V<:AbstractSimplex{0, T}}

A filtration is used to find the edges in filtration and to determine diameters of
simplices.

`T` is the distance type, accessible by `dist_type` and S is the edge type, accesible by
`vertex_type`.

# Interface

* [`n_vertices(::AbstractFiltration)`](@ref)
* [`edges(::AbstractFiltration)`](@ref)
* [`diam(::AbstractFiltration, vs)`](@ref)
* [`diam(::AbstractFiltration, ::AbstractSimplex, ::Any, ::Any)`](@ref)
* [`birth(::AbstractFiltration, v)`](@ref) - optional, defaults to returning `zero(T)`.
* [`threshold(::AbstractFiltration)`](@ref) - optional, defaults to returning `missing`.
"""
abstract type AbstractFiltration{T, V<:AbstractSimplex{0, T}} end

function Base.show(io::IO, flt::AbstractFiltration)
    print(io, typeof(flt), "(n_vertices=$(n_vertices(flt)))")
end

vertex_type(::AbstractFiltration{<:Any, V}) where V = V
edge_type(::AbstractFiltration{<:Any, V}) where V = coface_type(V)
dist_type(::AbstractFiltration{T}) where T = T

"""
    n_vertices(filtration::AbstractFiltration)

Return the number of vertices in `filtration`.
"""
n_vertices(::AbstractFiltration)

"""
    edges(filtration::AbstractFiltration)

Get edges in distance matrix in `filtration`, sorted by decresing length and increasing
combinatorial index. Edges should be of type `edge_type(filtration)`.
"""
edges(::AbstractFiltration)

"""
    diam(flt::AbstractFiltration, vertices)

Get the diameter of list of vertices i.e. diameter of simplex with `vertices`. Retrun
`missing` if `vertices` do not form a valid simplex.
"""
diam(::AbstractFiltration, ::Any)

"""
    diam(flt::AbstractFiltration, simplex, vertices, vertex)

Get the diameter of coface of a `Simplex` that is formed by adding `vertex` to `vertices`.
Return `missing` if new simplex is not valid.
"""
diam(::AbstractFiltration, ::AbstractSimplex, ::Any, ::Any)

"""
    birth(::AbstractFiltration, v)

Get the birth time of vertex `v`.
"""
birth(::AbstractFiltration{T}, _) where T = zero(T)

"""
    threshold(::AbstractFiltration)

Get the threshold of filtration. This is the maximum diameter a simplex in the filtration
can have. Used only for placing the infinity line in plotting.
"""
threshold(::AbstractFiltration) = missing
