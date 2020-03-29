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
    index(simplex::AbstractSimplex)

Get the combinatorial index of `simplex`. The index is equal to

```math
(i_d, i_{d-1}, ..., 1, 0) \\mapsto \\sum_{k=0^d} \\binom{i_k}{k + 1}.
```

    index(vertices, binomial)

Compute the index from a collection of `vertices`. Vertices must be in descending order.
"""
index(vertices, binomial) =
    sum(binomial(vertices[end - l + 1] - 1, l) for l in eachindex(vertices)) + 1

"""
    Simplex{M}

A simplex is represented by its combinatorial index (accessed by the `index` function) and
coefficient value (accessed by `coef`). The type parameter `M` represents the modulus of the
field of coefficients.

Note that the coefficient and value are stored in a single 8-byte word, so the range of
possible indices is slightly smaller than `typemax(Int64)`, depending on the number of bits
needed to represent `M`. A simplex has no information about its dimension.

# Constructor:

    Simplex{M}(index::Integer, value::Integer)
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
Simplex(x, coef, modulus) = Simplex{modulus}(x, coef)
# This constructor is for debugging and testing only.
Simplex{M}(vertices::AbstractVector, coef) where M =
    Simplex{M}(index(vertices, binomial), coef)

@generated function index(sx::Simplex{M}) where M
    bits = n_bits(M)
    :(reinterpret(Int64, sx) >> $bits)
end
@generated function coef(sx::Simplex{M}) where M
    bits = n_bits(M)
    :(reinterpret(Int64, sx) & (1 << $bits - 1))
end
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
"""
struct DiameterSimplex{M, T} <: AbstractSimplex{M}
    diam    ::T
    simplex ::Simplex{M}
end
DiameterSimplex{M}(diam, x, coef) where M =
    DiameterSimplex(diam, Simplex{M}(x, coef))
DiameterSimplex(diam, x, coef, modulus) =
    DiameterSimplex(diam, Simplex{modulus}(x, coef))

index(sx::DiameterSimplex) =
    index(sx.simplex)
coef(sx::DiameterSimplex) =
    coef(sx.simplex)
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
    get_vertices!(buffer, simplex::AbstractSimplex, dim, n_max, binomial)

Copy vertices of `dim`-dimensional `simplex` to `buffer`. `n_max` is the largest possible
vertex index i.e. the number of points in data set and `binomial` is a function or callable
object that is used to compute binomial coefficients.
"""
function get_vertices!(buff, sx::AbstractSimplex, dim, n_max, binomial)
    resize!(buff, dim + 1)
    idx = index(sx) - 1
    for (i, k) in enumerate(dim+1:-1:1)
        v = findlast(x -> binomial(x, k) ≤ idx, k - 1, n_max)
        buff[i] = v + 1
        idx -= binomial(v, k)
        n_max = v - 1
    end
    buff
end
