# prime fields =========================================================================== #
"""
    is_prime(n)

Return `true` if `n` is a prime number.
"""
function is_prime(n::Int)
    if iseven(n) || n < 2
        return n == 2
    else
        p = 3
        q = n ÷ p
        while p ≤ q
            iszero(n % p) && return false
            p += 2
            q = n ÷ p
        end
        return true
    end
end
# The following causes an error if you try to use anything other than Int as a modulus.
is_prime(::Any) = false

"""
    mod_prime(i, ::Val{M})

Like `mod`, but with prime `M`.
"""
function mod_prime(i, ::Val{M}) where {M}
    is_prime(M) || throw(DomainError(M, "modulus must be a prime number"))
    i = i % M
    return i + ifelse(signbit(i), M, 0)
end
mod_prime(i, ::Val{2}) = i & 1

"""
    Mod{M} <: Integer

`Mod{M}` is the default field used by Ripserer. It is a representation of a finite field
``\\mathbb{Z}_M``, integers modulo small, prime `M`. Supports field arithmetic and can be
converted to integer with `Int`.

Its values are not comparable on purpose.

# Example

```jldoctest
julia> Mod{3}(5)
2 mod 3

julia> Mod{3}(5) + 1
0 mod 3

```
"""
struct Mod{M} <: Integer
    value::Int

    # Check mod allows construction when you know you don't need to mod the number.
    function Mod{M}(value::Integer, check_mod=true) where {M}
        if check_mod
            return new{M}(mod_prime(value, Val(M)))
        else
            return new{M}(value)
        end
    end
end
Mod{M}(i::Mod{M}) where {M} = i

Base.Int(i::Mod) = i.value

Base.show(io::IO, i::Mod{M}) where {M} = print(io, Int(i), " mod ", M)

for op in (:+, :-, :*)
    @eval (Base.$op)(i::Mod{M}, j::Mod{M}) where {M} = Mod{M}($op(Int(i), Int(j)))
end

Base.:/(i::Mod{M}, j::Mod{M}) where {M} = i * inv(j)
Base.:-(i::Mod{M}) where {M} = Mod{M}(M - Int(i), false)
Base.zero(::Type{Mod{M}}) where {M} = Mod{M}(0, false)
Base.one(::Type{Mod{M}}) where {M} = Mod{M}(1, false)
Base.sign(i::M) where {M<:Mod} = ifelse(iszero(i), zero(M), one(M))

Base.promote_rule(::Type{Mod{M}}, ::Type{<:Mod}) where {M} = Union{}
function Base.promote_rule(::Type{Mod{M}}, ::Type{I}) where {M,I<:Integer}
    if Base.promote_type(I, Int128) === Int128 || Base.promote_type(I, Int128) === UInt128
        return Mod{M}
    else
        return Union{}
    end
end

Base.inv(i::Mod{M}) where {M} = Mod{M}(invmod(Int(i), M), false)
