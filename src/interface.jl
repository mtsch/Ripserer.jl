# abstract simplex ======================================================================= #
"""
    AbstractSimplex{D, C, T}

An abstract type for representing simplices. A simplex is represented by its dimension,
diameter, combinatorial index and coefficient value. It does not need to hold information
about its the vertices it includes, since they can be recomputed from the index and
dimension.

`D` is the dimension, `T` is the type of distance and `C` is the coefficient type. `D` is
accessible by `dim(::AbstractSimplex)`.

# Interface

* [`index(::AbstractSimplex)`](@ref)
* [`coef(::AbstractSimplex)`](@ref)
* [`set_coef(::AbstractSimplex, ::Any)`](@ref)
* [`diam(::AbstractSimplex)`](@ref)
* [`coface_type(::AbstractSimplex)`](@ref)
* [`vertices(::AbstractSimplex)`](@ref) - optional, comes with a default implementation.
* [`coboundary(::AbstractFiltration, ::AbstractSimplex)`](@ref) - optional, comes with a
  default implementation.
"""
abstract type AbstractSimplex{D, C, T} end

"""
    coef(simplex::AbstractSimplex)

Get the coefficient value of `simplex`.
"""
coef(::AbstractSimplex)

"""
    set_coef(simplex::AbstractSimplex, value)

Return new `simplex` of the same type with new coefficient `value`.
"""
set_coef(::AbstractSimplex, ::Any)

"""
    diam(simplex::AbstractSimplex)

Get the diameter of `simplex`.
"""
diam(::AbstractSimplex)

"""
    index(simplex::AbstractSimplex)

Get the combinatorial index of `simplex`.
"""
index(::AbstractSimplex)

"""
    coface_type(::AbstractSimplex)
    coface_type(::Type{<:AbstractSimplex})

Get the type of coface a simplex hax. For a `D`-dimensional simplex, this is usually its
`D+1`-dimensional counterpart. Only the `coface_type(::Type{<:AbstractSimplex})` method
needs to be implemented.
"""
coface_type(sx::AbstractSimplex) =
    coface_type(typeof(sx))

"""
    dim(::AbstractSimplex)
    dim(::Type{AbstractSimplex})

Get the dimension of simplex i.e. the value of `D`.
"""
dim(::Type{<:AbstractSimplex{D}}) where D =
    D
dim(sx::AbstractSimplex) =
    dim(typeof(sx))

# filtrations ============================================================================ #
"""
    AbstractFiltration{T, S<:AbstractSimplex{C, T}}

A filtration is used to find the edges in filtration and to determine diameters of
simplices.

`T` is the distance type, accessible by `dist_type` and S is the edge type, accesible by
`edge_type`.

# Interface

* [`n_vertices(::AbstractFiltration)`](@ref)
* [`edges(::AbstractFiltration)`](@ref)
* [`diam(::AbstractFiltration, vs)`](@ref)
* [`diam(::AbstractFiltration, ::AbstractSimplex, ::Any, ::Any)`](@ref)
* [`SparseArrays.issparse(::Type{A}) where A<:AbstractFiltration`](@ref) - optional defaults
  to `false`.
"""
abstract type AbstractFiltration{T, S<:AbstractSimplex{1, <:Any, T}} end

function Base.show(io::IO, flt::AbstractFiltration)
    print(io, typeof(flt), "(n_vertices=$(n_vertices(flt)))")
end

edge_type(::AbstractFiltration{<:Any, S}) where S =
    S
dist_type(::AbstractFiltration{T}) where T =
    T

"""
    n_vertices(filtration::AbstractFiltration)

Number of vertices in `filtration`.
"""
n_vertices(::AbstractFiltration)

"""
    SparseArrays.issparse(::Type{A}) where A<:AbstractFiltration

Return true if `A` is a sparse filtration. A filtration should be sparse if most simplices
are to be skipped. Defaults to `false`.
"""
SparseArrays.issparse(flt::AbstractFiltration) =
    issparse(typeof(flt))
SparseArrays.issparse(::Type{A}) where A<:AbstractFiltration =
    false

"""
    edges(filtration::AbstractFiltration)

Get edges in distance matrix in `filtration`, sorted by decresing length and increasing
combinatorial index. Edges should be of type `edge_type(filtration)`.
"""
edges(::AbstractFiltration)

"""
    diam(flt::AbstractFiltration, vertices)

Get the diameter of list of vertices i.e. diameter of simplex with `vertices`.
"""
diam(::AbstractFiltration, ::Any)
"""
    diam(flt::AbstractFiltration, simplex, vertices, vertex)

Get the diameter of coface of `simplex` that is formed by adding `vertex` to `vertices`.
"""
diam(::AbstractFiltration, ::AbstractSimplex, ::Any, ::Any)
