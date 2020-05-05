# prime fields =========================================================================== #
"""
    is_prime(n)

Return `true` if `n` is a prime number.
"""
@pure function is_prime(n::Int)
    if iseven(n) || n < 2
        n == 2
    else
        p = 3
        q = n / p
        while p ≤ q
            iszero(n % p) && return false
            p += 2
            q = n / p
        end
        true
    end
end
@pure is_prime(::Any) =
    false

"""
    mod_prime(i, ::Val{M})

Like `mod`, but with prime `M`.
"""
function mod_prime(i, ::Val{M}) where M
    is_prime(M) || throw(DomainError(M, "modulus must be a prime number"))
    i = i % M
    i + ifelse(signbit(i), M, 0)
end
mod_prime(i, ::Val{2}) =
    i & 1

"""
    PrimeField{M} <: Integer

Representation of finite field ``\\mathbb{Z}_M``, integers modulo `M`. Supports integer
arithmetic and can be converted to integer with `Int`.
"""
struct PrimeField{M} <: Integer
    value::Int

    function PrimeField{M}(value::Integer; check_mod=true) where M
        if check_mod
            new{M}(mod_prime(value, Val(M)))
        else
            new{M}(value)
        end
    end
end
PrimeField{M}(i::PrimeField{M}) where M =
    i

mod_prime(i::PrimeField{M}, ::Val{M}) where M =
    i
# prevent ambiguity.
mod_prime(i::PrimeField{2}, ::Val{2}) =
    i

Base.Int(i::PrimeField) =
    i.value

Base.show(io::IO, i::PrimeField{M}) where M =
    print(io, Int(i), " mod ", M)

for op in (:+, :-, :*)
    @eval (Base.$op)(i::PrimeField{M}, j::PrimeField{M}) where M =
        PrimeField{M}($op(Int(i), Int(j)))
end

Base.:/(i::PrimeField{M}, j::PrimeField{M}) where M =
    i * inv(j)
Base.:-(i::PrimeField{M}) where M =
    PrimeField{M}(M - Int(i), check_mod=false)
Base.zero(::Type{PrimeField{M}}) where M =
    PrimeField{M}(0, check_mod=false)
Base.one(::Type{PrimeField{M}}) where M =
    PrimeField{M}(1, check_mod=false)
Base.promote_rule(::Type{PrimeField{M}}, ::Type{Int}) where {M} =
    PrimeField{M}

# Idea: precompute inverses and generate a function with the inverses hard-coded.
@generated function Base.inv(i::PrimeField{M}) where M
    err_check = quote
        iszero(i) && throw(DivideError())
    end
    if M != 2
        inverse_arr = fill(PrimeField{M}(0), M-1)
        inverse_arr[1] = PrimeField{M}(1)
        for i in 2:M-1
            inverse_arr[i] = PrimeField{M}(invmod(i, M))
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

# default simplex ======================================================================== #
"""
    n_bits(M)

Get numer of bits needed to represent number mod `M`.
"""
@pure n_bits(M::Int) =
    floor(Int, log2(M-1)) + 1

"""
    Simplex{D, M, T} <: AbstractSimplex{D, PrimeField{M}, T}

The vanilla simplex type with coefficient values from ``\\mathbb{Z}_M``, integers modulo
`M`, where `M` is prime.

Because the algorithm needs to store and insane number of simplices, we pack `index` and
`coef` into a single value of type `U<:Unsigned`. The assumption here is that the modulus is
very small and `coef` only takes a few bits to store.

# Constructor

    Simplex{D, M}(::T, index::Integer, coef::Integer)

# Examples

```jldoctest
julia> Simplex{2, 2}(1, 2, 3)
2-dim Simplex{2}(1, 2, 1):
  (4, 2, 1)

julia> Simplex{10, 3}(1.0, Int128(10), 2)
4-dim Simplex{3}(1.0, 10, 2) with UInt128 index:
  (12, 11, 10, 9, 8, 7, 6, 5, 4, 2, 1)
`˙`
"""
struct Simplex{D, M, T, U<:Unsigned} <: AbstractSimplex{D, PrimeField{M}, T}
    diam       ::T
    index_coef ::U

    function Simplex{D, M, T, U}(
        diam::T, index::I, coef::Integer
    ) where {D, M, T, I<:Integer, U}
        is_prime(M) || throw(DomainError(M, "modulus must be a prime number"))
        D ≥ 0 || throw(DomainError(D, "dimension must be a non-negative integer"))
        U === unsigned(I) || throw(DomainError(U, "type parameters must match"))
        bits = n_bits(M)
        new{D, M, T, unsigned(I)}(
            diam, unsigned(index) << bits + I(mod_prime(coef, Val(M)))
        )
    end
end

Simplex{D, M}(diam::T, index_or_vertices, coef::Integer) where {D, M, T} =
    Simplex{D, M, T}(diam, index_or_vertices, coef)
Simplex{D, M, T}(diam, index::I, coef::Integer) where {D, M, T, I<:Integer} =
    Simplex{D, M, T, unsigned(I)}(T(diam), index, coef)
Simplex{D, M, T}(diam, vertices, coef::Integer) where {D, M, T} =
    Simplex{D, M, T}(T(diam), index(vertices), coef)

function Base.show(io::IO, ::MIME"text/plain", sx::Simplex{D, M, <:Any, U}) where {D, M, U}
    print(io, D, "-dim Simplex{", M, "}", (diam(sx), index(sx), Int(coef(sx))))
    if !(U <: UInt)
        print(io, " with ", U, " index")
    end
    print(io, ":\n  ", vertices(sx))
end

Base.show(io::IO, sx::Simplex{D, M}) where {D, M} =
    print(io, "Simplex{", D, ", ", M, "}", (diam(sx), vertices(sx), Int(coef(sx))))

# Interface implementation =============================================================== #
index(sx::Simplex{<:Any, M}) where M =
    signed(sx.index_coef >> n_bits(M))

function coef(sx::Simplex{<:Any, M}) where M
    mask = 1 << n_bits(M) - 1
    PrimeField{M}(sx.index_coef & mask, check_mod=false)
end

diam(sx::Simplex) =
    sx.diam

set_coef(sx::Simplex{D, M, T}, coef) where {D, M, T} =
    Simplex{D, M, T}(diam(sx), index(sx), coef)

@pure coface_type(::Type{Simplex{D, M, T, U}}) where {D, M, T, U} =
    Simplex{D+1, M, T, U}
