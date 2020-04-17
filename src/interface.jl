"""
    AbstractSimplex{C, T}

An abstract type for representing simplices. It is represented by a combinatorial index
and does not need to hold information about its dimension or the vertices it includes.
`T` is the type of distance and `C` is the coefficient type.

# Interface

    index(sx)::Int

    coef(sx)::C

    set_coef(sx)::typeof(sx)

    diam(sx)::T
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

    diam(flt::AbstractFiltration, vertices)

    diam(flt::AbstractFiltration, vertices, vertex)

Get the diameter of `simplex` or list of vertices. If additional `vertex` is given, only
calculate max distance from `vertices` to `vertex`.
"""
diam

"""
    index(simplex::AbstractSimplex)

Get the combinatorial index of `simplex`. The index is equal to

```math
(i_d, i_{d-1}, ..., 1) \\mapsto \\sum_{k=1}^{d+1} \\binom{i_k - 1}{k}.
```

    index(filtration::AbstractFiltration, vertices)

Compute the index from a collection of `vertices`. Vertices must be in descending order.
"""
index

"""
    AbstractFiltration{T, S<:AbstractSimplex{C, T}}

An abstract type that holds information about the distances between vertices and the simplex
type.

# Interface

    Base.length(::AbstractFiltration)::Int

    dist(::AbstractFiltration, ::Int, ::Int)::T

    edges(::AbstractFiltration)::iteratble of Tuple{T, {Int, Int}}

    Base.binomial(::AbstractFiltration, n, k)::Int (optional)

    dim_max(::AbstractFiltration)::Int

    threshold(::AbstractFiltration)::T
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

Get edges in distance matrix in `filtration`,
sorted by decresing length and increasing index.
"""
edges

Base.binomial(::AbstractFiltration, n, k) =
    binomial(n, k)

"""
    dim_max(flt::AbstractFiltration)

Get the maximum dimension of simplices in `flt`.
"""
dim_max

"""
    threshold(flt::AbstractFiltration)

Get the threshold of `flt`. Simplices with diameter strictly larger than this value will be
ignored.
"""
threshold(flt::AbstractFiltration) =
    infinity(flt)
