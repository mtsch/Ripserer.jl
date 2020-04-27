# abstract simplex ======================================================================= #
"""
    AbstractSimplex{D, C, T}

An abstract type for representing simplices. A simplex is represented by its dimension,
diameter, combinatorial index and coefficient value. It does not need to hold information
about its the vertices it includes, since they can be recomputed from the index and
dimension.

`D` is the dimension, `T` is the type of distance and `C` is the coefficient type.

# Interface

* [`index(::AbstractSimplex)`](@ref)
* [`coef(::AbstractSimplex)`](@ref)
* [`set_coef(::AbstractSimplex, ::Any)`](@ref)
* [`diam(::AbstractSimplex)`](@ref)
* [`vertices(::AbstractSimplex)`](@ref)
"""
abstract type AbstractSimplex{D, C, T} end

"""
    coef(simplex::AbstractSimplex)

Get the coefficient value of `simplex`.
"""
coef(::AbstractSimplex)

"""
    set_coef(simplex::AbstractSimplex, value)

Return new `simplex` with new coefficient `value`.
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

Get the type of coface a simplex hax. For a `D`-dimensional simplex, this is usually its
`D+1`-dimensional counterpart.
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

An abstract type that holds information about the distances between vertices and the simplex
type.

# Interface

* `n_vertices(::AbstractFiltration)` - return number of vertices in filtration.
* `edges(::AbstractFiltration)` - return all edges in filtration as `(l, (i, j))` where `l`
  is the edge length and `i` and `j` are its endpoints.
* `diam(::AbstractFiltration, vs)` - diameter of simplex with vertices in `vs`. Should
  return `Infinity()` if simplex is above threshold.
* `diam(::AbstractFiltration, sx::AbstractSimplex, vs, u)` - diameter of simplex `sx` with
  vertices in `vs` and an added vertex `u`. Should return `Infinity()` if simplex is above
  threshold.
* `SparseArrays.issparse(::Type{A}) where A<:AbstractFiltration` - optional, defaults to
  `false`. Should be `true` if most of the simplices are expected to be skipped.
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
diam(::AbstractFiltration, ::Any, ::Any, ::Any)
