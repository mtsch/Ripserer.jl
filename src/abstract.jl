"""
    AbstractSimplex{M, T}

An abstract type for representing simplices. It is represented by a combinatorial index
and does not need to hold information about its dimension or the vertices it includes.

# Interface

    index(sx)::Int

    coef(sx)::Int

    set_coef(sx)::typeof(sx)

    diam(sx)::T
"""
abstract type AbstractSimplex{M, T} end

"""
    coef(simplex::AbstractSimplex{M})

Get the coefficient value of `simplex`. The coefficient is always in the range of
`0 ≤ coef(simplex) < M`.
"""
coef

"""
    set_coef(simplex::AbstractSimplex, value)

Return new `simplex` with new coefficient `value`.
"""
set_coef

"""
    diam(sx::AbstractSimplex)

Get the diameter of simplex.
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
    AbstractFiltration{M, T, S<:AbstractSimplex{M, T}}

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
abstract type AbstractFiltration{M, T, S<:AbstractSimplex{M, T}} end

function Base.show(io::IO, flt::AbstractFiltration{M, T}) where {M, T}
    print(io, typeof(flt).name, "(length=$(length(flt)), modulus=$M")
    if threshold(flt) < typemax(T)
        print(io, ", threshold=$(threshold(flt))")
    end
    print(io, ", dim_max=$(dim_max(flt)), eltype=$(eltype(flt)))")
end

Base.eltype(::AbstractFiltration{M, T, S}) where {M, T, S} =
    S

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
threshold(flt::AbstractFiltration{M, T}) where {M, T} =
    typemax(T)

# implementation ========================================================================= #
Base.isless(sx1, sx2) =
    diam(sx1) < diam(sx2) || diam(sx1) == diam(sx2) && index(sx1) > index(sx2)

# simplex arithmetic ===================================================================== #
Base.:+(sx1::AbstractSimplex{M}, sx2::AbstractSimplex{M}) where M =
    set_coef(sx1, coef(sx1) + coef(sx2))
Base.:-(sx1::AbstractSimplex{M}, sx2::AbstractSimplex{M}) where M =
    set_coef(sx1, coef(sx1) - coef(sx2))
Base.:*(sx::AbstractSimplex, λ) =
    set_coef(sx, coef(sx) * λ)
Base.:*(λ, sx::AbstractSimplex) =
    set_coef(sx, λ * coef(sx))
Base.:-(sx::AbstractSimplex) =
    set_coef(sx, -coef(sx))
Base.:/(sx::AbstractSimplex{M}, λ) where M =
    set_coef(sx, coef(sx) * inv_mod(Val(M), λ))

"""
    inv_mod(::Val{M}, i)

Multiplicative inverse of `i` mod `M`.
"""
# Idea: precompute inverses and generate a function with the inverses hard-coded.
@generated function inv_mod(::Val{M}, i) where M
    err_check = quote
        i == 0 && throw(DivideError())
    end
    if M > 2
        isprime(M) || throw(DomainError(M, "modulus not prime"))

        inverse_arr = fill(0, M-1)
        inverse_arr[1] = 1
        for i in 2:M-1
            inverse_arr[i] = M - (inverse_arr[M % i] * floor(Int, M / i)) % M;
        end
        inverse = (inverse_arr...,)

        quote
            $err_check
            @inbounds $inverse[i]
        end
    else
        quote
            $err_check
            i
        end
    end
end

function index(flt::AbstractFiltration, vertices)
    res = 0
    for l in eachindex(vertices)
        res += binomial(flt, vertices[end - l + 1] - 1, l)
    end
    res + 1
end

function diam(flt::AbstractFiltration{M, T}, vertices) where {M, T}
    n = length(vertices)
    res = typemin(T)
    for i in 1:n, j in i+1:n
        d = dist(flt, vertices[j], vertices[i])
        if d == 0
            return typemax(T)
        else
            res = max(res, d)
        end
    end
    res
end

"""
    max_dist(filtration, vertices, vertex)

Get the maximum distance from `vertices` to `vertex`.
"""
function max_dist(flt::AbstractFiltration{M, T}, us, v::Integer) where {M, T}
    res = typemin(T)
    for u in us
        res = max(res, dist(flt, u, v))
    end
    res
end

"""
    find_max_vertex(filtration, idx, k)

Find largest vertex index of vertex for which `binomial(i, k) ≤ idx` holds.
"""
function find_max_vertex(flt::AbstractFiltration, idx, k)
    top = length(flt)
    bot = k - 1
    if !(binomial(flt, top, k) ≤ idx)
        count = top - bot
        while count > 0
            step = fld(count, 2)
            mid = top - step
            if !(binomial(flt, mid, k) ≤ idx)
                top = mid - 1
                count -= step + 1
            else
                count = step
            end
        end
    end
    top
end

"""
    get_vertices!(filtration, index)

Copy vertices of simplex with `index` to `filtration`'s vertex cache.
"""
function get_vertices!(flt::AbstractFiltration, index, dim)
    resize!(flt.vertex_cache, dim + 1)
    index = index - 1
    for (i, k) in enumerate(dim+1:-1:1)
        v = find_max_vertex(flt, index, k)

        flt.vertex_cache[i] = v + 1
        index -= binomial(flt, v, k)
        n_max = v - 1
    end
    flt.vertex_cache
end

"""
    vertices(filtration, simplex, dim)

    vertices(filtration, index, dim)

Get vertices of `simplex`. Vertices are only recomputed when the vertex cache in
`filtration` is invalid.
"""
vertices(flt::AbstractFiltration{M}, sx::AbstractSimplex{M}, dim) where M =
    vertices(flt, index(sx), dim)
function vertices(flt::AbstractFiltration, idx::Integer, dim)
    # Calculating index from vertices is so much faster that this is worth doing.
    if length(flt.vertex_cache) != dim+1 || index(flt, flt.vertex_cache) != idx
        get_vertices!(flt, idx, dim)
    end
    flt.vertex_cache
end

# coboundary ============================================================================= #
struct CoboundaryIterator{M, T, S<:AbstractSimplex{M, T}, F<:AbstractFiltration{M, T, S}}
    filtration  ::F
    simplex     ::S
    dim         ::Int
end

Base.IteratorSize(::Type{CoboundaryIterator}) =
    Base.SizeUnknown()
Base.IteratorEltype(::Type{CoboundaryIterator}) =
    Base.HasEltype()
Base.eltype(::Type{CoboundaryIterator{M, T, S}}) where {M, T, S} =
    S

"""
    coboundary(filtration, simplex, dim)

Return an iterator that iterates over all cofaces of `simplex` of dimension `dim + 1` in
decreasing order by index.
"""
coboundary(flt::AbstractFiltration, simplex::AbstractSimplex, dim) =
    CoboundaryIterator(flt, simplex, dim)

function Base.iterate(ci::CoboundaryIterator{M},
                      st = (length(ci.filtration),
                            ci.dim + 1,
                            index(ci.simplex) - 1,
                            0)) where M
    v, k, idx_below, idx_above = st
    v -= 1
    while v > 0 && v >= k && binomial(ci.filtration, v, k) <= idx_below
        idx_below -= binomial(ci.filtration, v, k)
        idx_above += binomial(ci.filtration, v, k + 1)
        v -= 1; k -= 1
    end
    if v < k
        nothing
    else
        vxs = vertices(ci.filtration, ci.simplex, ci.dim)
        diameter = max(diam(ci.simplex), max_dist(ci.filtration, vxs, v+1))

        coefficient = (k % 2 == 1 ? 1 : M - 1) * coef(ci.simplex) % M
        new_index = idx_above + binomial(ci.filtration, v, k + 1) + idx_below + 1
        eltype(ci.filtration)(diameter, new_index, coefficient), (v, k, idx_below, idx_above)
    end
end
