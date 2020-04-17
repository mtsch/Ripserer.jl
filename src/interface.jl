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

"""
    AbstractFiltration{T, S<:AbstractSimplex{C, T}}

An abstract type that holds information about the distances between vertices and the simplex
type.

# Interface

* `Base.length(::AbstractFiltration)`
* `dist(::AbstractFiltration, ::Integer, ::Integer)`
* `edges(::AbstractFiltration)`
* `dim_max(::AbstractFiltration)`
* `diam(::AbstractFiltration, iterable)` - optional, defaults to diameter of vertex set.
* `Base.binomial(::AbstractFiltration, n, k)` - optional, but recommended.
* `threshold(::AbstractFiltration)` - optional, defaults to `typemax(T)`.
"""
abstract type AbstractFiltration{T, S<:AbstractSimplex{<:Any, T}} end

function Base.show(io::IO, flt::AbstractFiltration)
    print(io, typeof(flt), "(length=$(length(flt))")
    if threshold(flt) < infinity(flt)
        print(io, ", threshold=$(threshold(flt))")
    end
    print(io, ", dim_max=$(dim_max(flt)))")
end

Base.eltype(::AbstractFiltration{<:Any, S}) where S =
    S
disttype(::AbstractFiltration{T}) where T =
    T
infinity(::AbstractFiltration{T}) where T =
    typemax(T)

"""
    length(filtration::AbstractFiltration)

Number of vertices in `filtration`.
"""
Base.length(::AbstractFiltration)

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
    dist(filtration::AbstractFiltration, i, j)

Get the distance between vertex `i` and vertex `j`.
"""
dist

"""
    edges(filtration::AbstractFiltration)

Get edges in distance matrix in `filtration`, sorted by decresing length and increasing
combinatorial index.
"""
edges

"""
    binomial(filtration::AbstractFiltration, n, k)

An abstract filtration may have binomial coefficients precomputed for better performance.
"""
Base.binomial(flt::AbstractFiltration, n, k) =
    binomial(n, k)

"""
    dim_max(filtration::AbstractFiltration)

Get the maximum dimension of simplices in `filtration`.
"""
dim_max

"""
    threshold(flt::AbstractFiltration)

Get the threshold of `flt`. Simplices with diameter strictly larger than this value will be
ignored.
"""
threshold(flt::AbstractFiltration) =
    infinity(flt)

"""
    diam(flt::AbstractFiltration, vertices)
    diam(flt::AbstractFiltration, vertices, vertex)

Get the diameter of list of vertices i.e. diameter of simplex with `vertices`. If
additional `vertex` is given, only calculate max distance from `vertices` to `vertex`.
"""
diam(::AbstractFiltration, ::Any)
