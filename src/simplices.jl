"""
    AbstractSimplex{M, T}

An abstract type for representing simplices. A simplex must support the following functions:

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
    diam(sx::Simplex)

Get the diameter of simplex.
"""
diam

"""
    index(simplex::AbstractSimplex)

Get the combinatorial index of `simplex`. The index is equal to

```math
(i_d, i_{d-1}, ..., 1, 0) \\mapsto \\sum_{k=0^d} \\binom{i_k}{k + 1}.
```

    index(reduction_state, vertices)

Compute the index from a collection of `vertices`. Vertices must be in descending order.
"""
function index(st::ReductionState, vertices)
    res = 0
    for l in eachindex(vertices)
        res += binomial(st, vertices[end - l + 1] - 1, l)
    end
    res + 1
end

"""
    n_bits(::Val{M})

Get numer of bits needed to represent number mod `M`.
"""
n_bits(::Val{M}) where M =
    floor(Int, log2(M-1)) + 1

struct Simplex{M, T} <: AbstractSimplex{M, T}
    diam       ::T
    index_coef ::Int64
end

@generated function Simplex{M}(diam::T, index, coef) where {M, T}
    isprime(M) || throw(DomainError(M, "modulus not prime"))
    bits = n_bits(Val(M))
    :(Simplex{M, T}(diam, Int64(index) << $bits + mod(coef, $M)))
end
Simplex{M}(st::ReductionState{M}, diam, vertices, coef) where M =
    Simplex{M}(diam, index(st, vertices), coef)

index(sx::Simplex{M}) where M =
    sx.index_coef >> n_bits(Val(M))

coef(sx::Simplex{M}) where M =
    sx.index_coef & (1 << n_bits(Val(M)) - 1)

diam(sx::Simplex) =
    sx.diam

set_coef(sx::Simplex{M}, value) where M =
    Simplex{M}(diam(sx), index(sx), value)

Base.show(io::IO, sx::Simplex{M}) where M =
    print(io, "Simplex{", M, "}", (diam(sx), index(sx), coef(sx)));

"""
    SimplexComparer

Ordering on `AbstractSimplex` by

* increasing diameter,
* decreasing combinatorial index.
"""
struct SimplexComparer end

(::SimplexComparer)(sx1, sx2) =
    diam(sx1) < diam(sx2) || diam(sx1) == diam(sx2) && index(sx1) > index(sx2)

DataStructures.compare(sc::SimplexComparer, sx1, sx2) =
    sc(sx1, sx2)

function find_max_vertex(st, idx, k)
    top = n_vertices(st)
    bot = k - 1
    if !(binomial(st, top, k) ≤ idx)
        count = top - bot
        while count > 0
            step = fld(count, 2)
            mid = top - step
            if !(binomial(st, mid, k) ≤ idx)
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
    get_vertices!(reduction_state, simplex::AbstractSimple)

Copy vertices of `simplex` to `reduction_state`'s vertex cache.
"""
function get_vertices!(st::ReductionState, sx::AbstractSimplex, dim)
    resize!(st.vertex_cache, dim + 1)
    idx = index(sx) - 1
    for (i, k) in enumerate(dim+1:-1:1)
        v = find_max_vertex(st, idx, k)

        st.vertex_cache[i] = v + 1
        idx -= binomial(st, v, k)
        n_max = v - 1
    end
    st.vertex_cache
end

"""
    vertices(reduction_state, simplex)

Get vertices of `simplex`. Vertices are only recomputed when the vertex cache in
`reduction_state` is invalid.
"""
function vertices(st::ReductionState{M}, sx::AbstractSimplex{M}, dim) where M
    # Calculating index from vertices is so much faster that this is worth doing.
    if length(st.vertex_cache) != dim+1 || index(st, st.vertex_cache) != index(sx)
        get_vertices!(st, sx, dim)
    end
    st.vertex_cache
end

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
            $inverse[i]
        end
    else
        quote
            $err_check
            i
        end
    end
end
