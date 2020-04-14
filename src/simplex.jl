# prime fields =========================================================================== #
"""
    isprime(n)

Return `true` if `n` is a prime number.
"""
function isprime(n)
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

"""
    modprime(i, M)

Like `mod`, but with prime, i.e. positive `M`.
"""
function modprime(i, M)
    i = i % M
    i + ifelse(signbit(i), M, 0)
end

"""
    PrimeField{M} <: Integer

Representation of finite field, integers modulo `M`.
"""
primitive type PrimeField{M} <: Integer 64 end

@generated function PrimeField{M}(value::Integer) where M
    isprime(M) || throw(DomainError(M, "modulus not prime"))
    if M == 2
        :(reinterpret(PrimeField{$M}, value & 1))
    else
        :(reinterpret(PrimeField{$M}, modprime(value, $M)))
    end
end
PrimeField{M}(i::PrimeField{M}) where M =
    i

Base.show(io::IO, i::PrimeField{M}) where M =
    print(io, Int(i), " mod ", M)

Base.Int(i::PrimeField) =
    reinterpret(Int, i)

for op in (:+, :-, :*)
    @eval (Base.$op)(i::PrimeField{M}, j::PrimeField{M}) where M =
        PrimeField{M}($op(Int(i), Int(j)))
end

Base.:/(i::PrimeField{M}, j::PrimeField{M}) where M =
    i * inv(j)
Base.:-(i::PrimeField{M}) where M =
    reinterpret(PrimeField{M}, M - Int(i))
Base.zero(::Type{PrimeField{M}}) where M =
    reinterpret(PrimeField{M}, 0)
Base.one(::Type{PrimeField{M}}) where M =
    reinterpret(PrimeField{M}, 1)
Base.promote_rule(::Type{PrimeField{M}}, ::Type{Int}) where {M} =
    PrimeField{M}

# Idea: precompute inverses and generate a function with the inverses hard-coded.
@generated function Base.inv(i::PrimeField{M}) where M
    err_check = quote
        iszero(i) && throw(DivideError())
    end
    if M != 2
        isprime(M) || throw(DomainError(M, "modulus not prime"))

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

# simplex ================================================================================ #
"""
    n_bits(M)

Get numer of bits needed to represent number mod `M`.
"""
n_bits(M) =
    floor(Int, log2(M-1)) + 1

"""
    Simplex{M, T} <: AbstractSimplex{PrimeField{M}, T}

The vanilla simplex type with coefficient values from `Z_M`, integers modulo `M`.
`index` and `coef` are packed into a single `UInt64`.

# Constructor

    Simplex{M}(::T, index::Integer, coef::Integer)
"""
struct Simplex{M, T} <: AbstractSimplex{PrimeField{M}, T}
    diam       ::T
    index_coef ::UInt64

    @generated function Simplex{M, T}(diam::T, index::Integer, coef::Integer) where {M, T}
        isprime(M) || throw(DomainError(M, "modulus not prime"))
        bits = n_bits(M)
        Expr(:new, :(Simplex{$M, $T}),
             :diam,
             :(UInt64(index) << $bits + modprime(coef, $M)))
    end
    @generated function Simplex{M, T}(diam::T, index::Integer, coef::PrimeField{M}) where {M, T}
        bits = n_bits(M)
        Expr(:new, :(Simplex{$M, $T}),
             :diam,
             :(UInt64(index) << $bits + Int(coef)))
    end
end

Simplex{M}(diam::T, index, coef) where {M, T} =
    Simplex{M, T}(diam, index, coef)
Simplex{M}(flt::AbstractFiltration{T}, diam, vertices, coef) where {M, T} =
    Simplex{M, T}(diam, index(flt, vertices), coef)
Simplex{M, T}(flt::AbstractFiltration{T}, diam, vertices, coef) where {M, T} =
    Simplex{M, T}(diam, index(flt, vertices), coef)

@generated function index(sx::Simplex{M}) where M
    shift = n_bits(M)
    :(reinterpret(Int64, sx.index_coef >> $shift))
end

@generated function coef(sx::Simplex{M}) where M
    mask = 1 << n_bits(M) - 1
    :(reinterpret(PrimeField{$M}, sx.index_coef & $mask))
end

diam(sx::Simplex) =
    sx.diam

set_coef(sx::Simplex{M, T}, coef) where {M, T} =
    Simplex{M, T}(diam(sx), index(sx), coef)

Base.show(io::IO, sx::Simplex{M}) where M =
    print(io, "Simplex{", M, "}", (diam(sx), index(sx), Int(coef(sx))))