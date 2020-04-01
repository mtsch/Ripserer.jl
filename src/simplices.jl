"""
    AbstractSimplex{M}

An abstract type for representing simplices. A simplex must support the following functions:

    index(sx)::Int
    coef(sx)::Int
"""
abstract type AbstractSimplex{M} end

"""
    coef(simplex::AbstractSimplex{M})

Get the coefficient value of `simplex`. The coefficient is always in the range of
`0 ≤ coef(simplex) < M`.
"""
coef

"""
    set_coef(simplex::AbstractSimplex, val)

Return new `simplex` with new coefficietn value.
"""
set_coef

"""
    index(simplex::AbstractSimplex)

Get the combinatorial index of `simplex`. The index is equal to

```math
(i_d, i_{d-1}, ..., 1, 0) \\mapsto \\sum_{k=0^d} \\binom{i_k}{k + 1}.
```

    index(reduction_state, vertices)

Compute the index from a collection of `vertices`. Vertices must be in descending order.
"""
index(st::ReductionState, vertices) =
    sum(binomial(vertices[end - l + 1] - 1, l) for l in eachindex(vertices)) + 1

"""
    Simplex{M}

A simplex is represented by its combinatorial index (accessed by the `index` function) and
coefficient value (accessed by `coef`). The type parameter `M` represents the modulus of the
field of coefficients.

Note that the coefficient and value are stored in a single 8-byte word, so the range of
possible indices is slightly smaller than `typemax(Int64)`, depending on the number of bits
needed to represent `M`. A simplex has no information about its dimension.

# Constructors:

    Simplex{M}(index::Integer, value::Integer)

    Simplex{M}(st::ReductionState, vertices, value::Integer)
"""
primitive type Simplex{M} <: AbstractSimplex{M} 64 end

"""
    n_bits(M)

Return the number of bits needed to represent an integer mod `M`.
"""
n_bits(M) =
    floor(Int, log2(M-1)) + 1

@generated function Simplex{M}(index, coef) where M
    isprime(M) || throw(DomainError(M, "modulus not prime"))
    bits = n_bits(M)
    :(reinterpret(Simplex{M}, Int64(index) << $bits + mod(coef, $M)))
end
Simplex{M}(st::ReductionState, vertices, coef) where M =
    Simplex{M}(index(st, vertices), coef)

@generated function index(sx::Simplex{M}) where M
    bits = n_bits(M)
    :(reinterpret(Int64, sx) >> $bits)
end
@generated function coef(sx::Simplex{M}) where M
    bits = n_bits(M)
    :(reinterpret(Int64, sx) & (1 << $bits - 1))
end
set_coef(sx::Simplex{M}, value) where M =
    Simplex{M}(index(sx), value)
Base.show(io::IO, ent::Simplex{M}) where M =
    print(io, "Simplex{$M}$((index(ent), coef(ent)))");

"""
    DiameterSimplex{M, T}

A simplex with a diameter is represented by its diameter (accesse by the `diam` function),
combinatorial index (accessed by `index`) and coefficient value (accessed by `coef`). The
type parameter `M` represents the modulus of the field of coefficients and `T` represents
the type of diameter.

Note that the coefficient and value are stored in a single 8-byte word, so the range of
possible indices is slightly smaller than `typemax(Int64)`, depending on the number of bits
needed to represent `M`. A simplex has no information about its dimension.

# Constructor:

    DiameterSimplex{M}(diameter::T, index::Integer, value::Integer)

    DiameterSimplex{M}(st::ReductionState, diameter::T, vertices, value::Integer)
"""
struct DiameterSimplex{M, T} <: AbstractSimplex{M}
    diam    ::T
    simplex ::Simplex{M}
end
DiameterSimplex{M}(diam, x, coef) where M =
    DiameterSimplex(diam, Simplex{M}(x, coef))
DiameterSimplex{M}(st::ReductionState, diam, vertices, coef) where M =
    DiameterSimplex(diam, Simplex{M}(st, vertices, coef))

index(sx::DiameterSimplex) =
    index(sx.simplex)
coef(sx::DiameterSimplex) =
    coef(sx.simplex)
set_coef(sx::DiameterSimplex{M}, value) where M =
    DiameterSimplex{M}(diam(sx), index(sx), value)
"""
    diam(sx::DiameterSimplex)

Get the diameter of simplex.
"""
diam(sx::DiameterSimplex) =
    sx.diam
Base.show(io::IO, sx::DiameterSimplex) =
    print(io, "DiameterSimplex$((diam(sx), index(sx), coef(sx)))")

"""
    DiameterSimplexComparer

Ordering on `DiameterSimplex` by

* increasing diameter,
* decreasing combinatorial index.
"""
struct DiameterSimplexComparer end

(::DiameterSimplexComparer)(sx1, sx2) =
    diam(sx1) < diam(sx2) || diam(sx1) == diam(sx2) && index(sx1) > index(sx2)

DataStructures.compare(dsc::DiameterSimplexComparer, sx1, sx2) =
    dsc(sx1, sx2)

# Find largest integer i between bot and top, for which f(i) is true.
function Base.findlast(f, bot::Int, top::Int)
    if !f(top)
        count = top - bot
        while count > 0
            step = count ÷ 2
            mid = top - step
            if !f(mid)
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
function get_vertices!(st::ReductionState, sx::AbstractSimplex)
    d = dim(st)
    resize!(st.vertex_cache, d + 1)
    idx = index(sx) - 1
    for (i, k) in enumerate(d+1:-1:1)
        v = findlast(x -> binomial(st, x, k) ≤ idx, k - 1, n_vertices(st))
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
function vertices(st::ReductionState{M}, sx::AbstractSimplex{M}) where M
    # Calculating index from vertices is so much faster that this is worth doing.
    if length(st.vertex_cache) != dim(st)+1 || index(st, st.vertex_cache) != index(sx)
        get_vertices!(st, sx)
    end
    st.vertex_cache
end

# simplex arithmetic ===================================================================== #

# Note: mod is handled in set_coef.
for op in (:+, :-, :*)
    @eval function (Base.$op)(sx1::AbstractSimplex{M}, sx2::AbstractSimplex{M}) where M
        @boundscheck begin
            # index(sx1) == index(sx2) || throw(ArgumentError("simplex indices don't match"))
        end
        set_coef(sx1, $op(coef(sx1), coef(sx2)))
    end
end

Base.:/(sx1::AbstractSimplex{M}, sx2::AbstractSimplex{M}) where M =
    sx1 * inv(sx2)
Base.:-(sx::AbstractSimplex{M}) where M =
    set_coef(sx, -coef(sx))

# Idea: precompute inverses and generate a function with the inverses hard-coded.
@generated function Base.inv(sx::AbstractSimplex{M}) where M
    err_check = quote
        coef(sx) == 0 && throw(DomainError(0))
    end
    if M > 2
        inverse_arr = fill(0, M-1)
        inverse_arr[1] = 1
        for i in 2:M-1
            inverse_arr[i] = M - (inverse_arr[M % i] * floor(Int, M / i)) % M;
        end
        inverse = (inverse_arr...,)

        quote
            $err_check
            set_coef(sx, $inverse[coef(sx)])
        end
    else
        quote
            $err_check
            set_coef(sx, coef(sx))
        end
    end
end
