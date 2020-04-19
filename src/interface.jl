# simplices ============================================================================== #
"""
    AbstractSimplex{C, T}

An abstract type for representing simplices. A simplex is represented by its diameter,
combinatorial index and coefficient value. It does not need to hold information about its
dimension or the vertices it includes.

`T` is the type of distance and `C` is the coefficient type.

# Interface

* `index(::AbstractSimplex)`
* `coef(::AbstractSimplex)`
* `set_coef(::AbstractSimplex{C}, ::C)`
* `diam(::AbstractSimplex)`
"""
abstract type AbstractSimplex{C, T} end

"""
    coef(simplex::AbstractSimplex)

Get the coefficient value of `simplex`.
"""
coef

"""
    set_coef(simplex::AbstractSimplex, value)

Return new `simplex` with new coefficient `value`.
"""
set_coef

"""
    diam(simplex::AbstractSimplex)

Get the diameter of `simplex`.
"""
diam(::AbstractSimplex)

"""
    index(simplex::AbstractSimplex)

Get the combinatorial index of `simplex`. The index of is equal to

```math
(i_d, i_{d-1}, ..., 1) \\mapsto \\sum_{k=1}^{d+1} \\binom{i_k - 1}{k},
```
where `i_k` are the simplex vertex indices.
"""
index(::AbstractSimplex)

# simplex arithmetic ===================================================================== #
Base.isless(sx1, sx2) =
    diam(sx1) < diam(sx2) || diam(sx1) == diam(sx2) && index(sx1) > index(sx2)

Base.:+(sx1::A, sx2::A) where A<:AbstractSimplex =
    set_coef(sx1, coef(sx1) + coef(sx2))
Base.:-(sx1::A, sx2::A) where A<:AbstractSimplex =
    set_coef(sx1, coef(sx1) - coef(sx2))
Base.:*(sx::AbstractSimplex, λ::Number) =
    set_coef(sx, coef(sx) * λ)
Base.:*(λ::Number, sx::AbstractSimplex) =
    set_coef(sx, λ::Number * coef(sx))
Base.:-(sx::AbstractSimplex) =
    set_coef(sx, -coef(sx))
Base.:/(sx::AbstractSimplex{C}, λ::Number) where C =
    set_coef(sx, coef(sx) * inv(C(λ)))

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
abstract type AbstractFiltration{T, S<:AbstractSimplex{<:Any, T}} end

function Base.show(io::IO, flt::AbstractFiltration)
    print(io, typeof(flt), "(n_vertices=$(n_vertices(flt))")
    if threshold(flt) < ∞
        print(io, ", threshold=$(threshold(flt))")
    end
    println(io, ")")
end

Base.eltype(::AbstractFiltration{<:Any, S}) where S =
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
combinatorial index.
"""
edges

"""
    diam(flt::AbstractFiltration, vertices)
    diam(flt::AbstractFiltration, vertices, vertex)

Get the diameter of list of vertices i.e. diameter of simplex with `vertices`. If
additional `vertex` is given, only calculate max distance from `vertices` to `vertex`.
"""
diam(::AbstractFiltration, ::Any)
