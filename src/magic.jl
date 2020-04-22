# Beware: black magic ahead.
using Ripserer
import Ripserer: index, vertices, coboundary

# prime fields =========================================================================== #
"""
    is_prime(n)

Return `true` if `n` is a prime number.
"""
function is_prime(n)
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
    prime_check(::Val{M})

Check that type parameter `M` is prime at compile time.
"""
@generated function prime_check(::Val{M}) where M
    is_prime(M) || throw(DomainError(M, "modulus not prime"))
    nothing
end

"""
    mod_prime(i, ::Val{M})

Like `mod`, but with prime `M`.
"""
function mod_prime(i, ::Val{M}) where M
    prime_check(Val(M))
    i = i % M
    i + ifelse(signbit(i), M, 0)
end
mod_prime(i, ::Val{2}) =
    i & 1

"""
    PrimeField{M} <: Integer

Representation of finite field, integers modulo `M`.
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
Base.promote_rule(::Type{PrimeField{M}}, ::Type{Int}) where {D, M} =
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

# magic simplex ========================================================================== #
"""
    n_bits(M)

Get numer of bits needed to represent number mod `M`.
"""
n_bits(M) =
    floor(Int, log2(M-1)) + 1

"""
    MagicSimplex{D, M, T} <: AbstractSimplex{PrimeField{M}, T}

The vanilla simplex type with coefficient values from ``\\mathbb{Z}_M``, integers modulo
prime `M`. `index` and `coef` are packed into a single `UInt64`.

# Constructor

    MagicSimplex{D, M}(::T, index::Integer, coef::Integer)
"""
struct MagicSimplex{D, M, T} <: AbstractSimplex{PrimeField{M}, T}
    diam       ::T
    index_coef ::UInt64

    @generated function MagicSimplex{D, M, T}(
        diam::T, index::Integer, coef::Integer
    ) where {D, M, T}

        prime_check(Val(M))
        bits = n_bits(M)
        Expr(:new, :(MagicSimplex{D, M, T}),
             :diam,
             :(UInt64(index) << $bits + mod_prime(coef, Val(M))))
    end
    @generated function MagicSimplex{D, M, T}(
        diam::T, index::Integer, coef::PrimeField{M}
    ) where {D, M, T}

        bits = n_bits(M)
        Expr(:new, :(MagicSimplex{D, M, T}),
             :diam,
             :(UInt64(index) << $bits + Int(coef)))
    end
end

MagicSimplex{D, M}(diam::T, index, coef) where {D, M, T} =
    MagicSimplex{D, M, T}(diam, index, coef)

@generated function index(sx::MagicSimplex{<:Any, M}) where M
    shift = n_bits(M)
    :(reinterpret(Int64, sx.index_coef >> $shift))
end

@generated function coef(sx::MagicSimplex{<:Any, M}) where M
    mask = 1 << n_bits(M) - 1
    :(PrimeField{$M}(sx.index_coef & $mask, check_mod=false))
end

diam(sx::MagicSimplex) =
    sx.diam

set_coef(sx::MagicSimplex{D, M, T}, coef) where {D, M, T} =
    MagicSimplex{D, M, T}(diam(sx), index(sx), coef)

Base.show(io::IO, sx::MagicSimplex{D, M}) where {D, M} =
    print(io, "MagicSimplex{", D, ", ", M, "}", (diam(sx), index(sx), Int(coef(sx))))

# magical stuff ========================================================================== #
magic_binomial(_, ::Val{0}) = 1
magic_binomial(n, ::Val{1}) = n

@generated function magic_binomial(n, ::Val{k}) where k
    coefficients = [0]
    k_minus_1 = k - 1

    quote
        n_coef = length($coefficients)
        if n_coef < n + 1
            resize!($coefficients, n)
            for i in n_coef+1:n
                $coefficients[i] = magic_binomial(i-1, Val($k_minus_1)) + $coefficients[i-1]
            end
        end
        @inbounds $coefficients[n]
    end
end

# Generate code of the form
# vk   = find_max_vertex(index, Val(k))
# vk-1 = find_max_vertex(index, Val(k-1), vk)
# ...
# v1 = find_max_vertex(index, Val(3), v2)
# v0 = find_max_vertex(index, Val(3), v1)
# (vk, ..., v0)
@generated function vertices(index, ::Val{dim}) where dim
    vars = Symbol[Symbol("v", k) for k in dim:-1:0]
    expr = quote
        index -= 1
        $(vars[1]) = find_max_vertex(index, Val($dim+1))
        index -= magic_binomial($(vars[1]), Val($dim+1))
    end

    for (i, k) in enumerate(dim:-1:1)
        expr = quote
            $expr
            $(vars[i+1]) = find_max_vertex(index, Val($k), $(vars[i]) - 1)
            index -= magic_binomial($(vars[i+1]), Val($k))
        end
    end
    quote
        $expr
        tuple($(vars...)) .+ 1
    end
end

vertices(sx::MagicSimplex{D}) where D =
    vertices(index(sx), Val(D))

function find_max_vertex(idx, ::Val{k}) where k
    n_max = 10
    while magic_binomial(n_max, Val(k)) ≤ idx
        n_max <<= 1
    end
    find_max_vertex(idx, Val(k), n_max)
end

function find_max_vertex(idx, ::Val{k}, n_max) where k
    hi = n_max + 1
    lo = k - 1
    while lo < hi - 1
        m = lo + ((hi - lo) >>> 0x01)
        if magic_binomial(m, Val(k)) ≤ idx
            lo = m
        else
            hi = m
        end
    end
    lo
end

# generate code of the form
# 1 + magic_binomial(vertices[1] - 1, Val(k))
#   + magic_binomial(vertices[2] - 1, Val(k-1))
#   ...
#   + magic_binomial(vertices[k] - 1, Val(1))
@generated function index(vertices::NTuple{k}) where k
    expr = quote
        1 + magic_binomial(vertices[1]-1, Val($k))
    end
    for i in 2:k
        expr = quote
            $expr + magic_binomial(vertices[$i]-1, Val($(k - i + 1)))
        end
    end
    expr
end

struct Coboundary{D1, D2, F, S<:MagicSimplex{D1}}
    filtration ::F
    simplex    ::S
    vertices   ::NTuple{D2, Int}
end

coboundary(filtration, simplex, dim) =
    Coboundary{dim, dim+1, typeof(filtration), typeof(simplex)}(
        filtration,
        simplex,
        vertices(simplex)
    )

function Base.iterate(ci::Coboundary{D1}, (v, k)=(n_vertices(filtration), D1+1)) where D1
    diameter = ∞
    @inbounds while diameter == ∞
        while v > 1 && v in ci.vertices
            v -= 1
            k -= 1
        end
        diameter = diam(ci.filtration, ci.simplex, ci.vertices, v)
    end
    if diameter != ∞
        @assert k ≥ 0
        coefficient = ifelse(k % 2 == 1, -coef(ci.simplex), coef(ci.simplex))
        # todo: be smarter than sort...
        new_index = index(TupleTools.sort(tuple(ci.vertices..., v), rev=true))

        eltype(ci.filtration)(diameter, new_index, coefficient), v - 1
    else
        nothing
    end
end

# benching
