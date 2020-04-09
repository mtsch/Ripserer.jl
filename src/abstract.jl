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
    coef(simplex::AbstractSimplex)

Get the coefficient value of `simplex`. The coefficient is always in the range of
`0 ≤ coef(simplex) < M`.
"""
coef

"""
    set_coef(simplex::AbstractSimplex, val)

Return new `simplex` with new coefficient value.
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
(i_d, i_{d-1}, ..., 1, 0) \\mapsto \\sum_{k=0^d} \\binom{i_k}{k + 1}.
```

    index(complex::SimplicialComplex, vertices)

Compute the index from a collection of `vertices`. Vertices must be in descending order.
"""
index

"""
    SimplicialComplex{M, T, S<:AbstractSimplex{M, T}}

An abstract type that holds information about the distances between vertices and the simplex
type.

# Interface

    Base.length(::SimplicialComplex)::Int

    dist(::SimplicialComplex, ::Int, ::Int)::T

    edges(::SimplicialComplex)::iteratble of Tuple{T, {Int, Int}}

    Base.binomial(::SimplicialComplex, n, k)::Int (optional)

    dim_max(::SimplicialComplex)::Int
"""
abstract type SimplicialComplex{M, T, S<:AbstractSimplex{M, T}} end

Base.eltype(::SimplicialComplex{M, T, S}) where {M, T, S} = S

"""
    length(complex::SimplicialComplex)

Get number of vertices in `complex`.
"""
length

"""
    dist(complex::SimplicialComplex, i, j)

Get the distance between vertex `i` and vertex `j`.
"""
dist

"""
    edges(complex::SimplicialComplex)

Get edges in distance matrix in `complex`,
sorted by decresing length and increasing index.
"""
edges

Base.binomial(::SimplicialComplex, n, k) =
    binomial(n, k)

"""
    dim_max(scx::SimplicialComplex)

Get the maximum dimension of simplices in `scx`.
"""
dim_max(scx::SimplicialComplex) =
    scx.dim_max

"""
    threshold(scx::SimplicialComplex)

Get the threshold of `scx`. Simplices with diameter strictly larger than this value will be
ignored.
"""
threshold(scx::SimplicialComplex{M, T}) where {M, T} =
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

function index(scx::SimplicialComplex, vertices)
    res = 0
    for l in eachindex(vertices)
        res += binomial(scx, vertices[end - l + 1] - 1, l)
    end
    res + 1
end

function diam(scx::SimplicialComplex{M, T}, vertices) where {M, T}
    n = length(vertices)
    res = typemin(T)
    for i in 1:n, j in i+1:n
        d = dist(scx, vertices[j], vertices[i])
        if d == 0
            return typemax(T)
        else
            res = max(res, d)
        end
    end
    res
end

"""
    max_dist(complex, vertices, vertex)

Get the maximum distance from `vertices` to `vertex`.
"""
function max_dist(scx::SimplicialComplex{M, T}, us, v::Integer) where {M, T}
    res = typemin(T)
    for u in us
        res = max(res, dist(scx, u, v))
    end
    res
end

"""
    is_connected(complex, vertices, vertex)

Check if `vertex` is connected to `vertices` i.e. if the distance to all other vertices
is ≥ 0. If `vertex in vertices`, this function returns `false`.
"""
function is_connected(scx::SimplicialComplex, us, v)
    for u in us
        if iszero(dist(scx, u, v))
            return false
        end
    end
    true
end

function find_max_vertex(scx::SimplicialComplex, idx, k)
    top = length(scx)
    bot = k - 1
    if !(binomial(scx, top, k) ≤ idx)
        count = top - bot
        while count > 0
            step = fld(count, 2)
            mid = top - step
            if !(binomial(scx, mid, k) ≤ idx)
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
    get_vertices!(complex, index)

Copy vertices of simplex with `index` to `complex`'s vertex cache.
"""
function get_vertices!(scx::SimplicialComplex, index, dim)
    resize!(scx.vertex_cache, dim + 1)
    index = index - 1
    for (i, k) in enumerate(dim+1:-1:1)
        v = find_max_vertex(scx, index, k)

        scx.vertex_cache[i] = v + 1
        index -= binomial(scx, v, k)
        n_max = v - 1
    end
    scx.vertex_cache
end

"""
    vertices(complex, simplex, dim)

    vertices(complex, index, dim)

Get vertices of `simplex`. Vertices are only recomputed when the vertex cache in
`complex` is invalid.
"""
vertices(scx::SimplicialComplex{M}, sx::AbstractSimplex{M}, dim) where M =
    vertices(scx, index(sx), dim)
function vertices(scx::SimplicialComplex, idx::Integer, dim)
    # Calculating index from vertices is so much faster that this is worth doing.
    if length(scx.vertex_cache) != dim+1 || index(scx, scx.vertex_cache) != idx
        get_vertices!(scx, idx, dim)
    end
    scx.vertex_cache
end

# coboundary ============================================================================= #
struct CoboundaryIterator{M, T, S<:AbstractSimplex{M, T}, C<:SimplicialComplex{M, T, S}}
    complex     ::C
    simplex     ::S
    dim         ::Int
end

Base.IteratorSize(::Type{CoboundaryIterator}) =
    Base.SizeUnknown()
Base.IteratorEltype(::Type{CoboundaryIterator}) =
    Base.HasEltype()
Base.eltype(::Type{CoboundaryIterator{M, T, S}}) where {M, T, S} =
    S

coboundary(scx::SimplicialComplex, simplex::AbstractSimplex, dim) =
    CoboundaryIterator(scx, simplex, dim)

function Base.iterate(ci::CoboundaryIterator{M},
                      st = (length(ci.complex),
                            ci.dim + 1,
                            index(ci.simplex) - 1,
                            0)) where M
    v, k, idx_below, idx_above = st
    v -= 1
    while v > 0 && v >= k && binomial(ci.complex, v, k) <= idx_below
        idx_below -= binomial(ci.complex, v, k)
        idx_above += binomial(ci.complex, v, k + 1)
        v -= 1; k -= 1
    end
    if v < k
        nothing
    else
        vxs = vertices(ci.complex, ci.simplex, ci.dim)
        diameter = max(diam(ci.simplex), max_dist(ci.complex, vxs, v+1))

        coefficient = (k % 2 == 1 ? 1 : M - 1) * coef(ci.simplex) % M
        new_index = idx_above + binomial(ci.complex, v, k + 1) + idx_below + 1
        eltype(ci.complex)(diameter, new_index, coefficient), (v, k, idx_below, idx_above)
    end
end
