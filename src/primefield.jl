# prime fields =========================================================================== #
"""
    is_prime(n)

Return `true` if `n` is a prime number.
"""
@pure function is_prime(n::Int)
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
@pure is_prime(::Any) = false

"""
    mod_prime(i, ::Val{M})

Like `mod`, but with prime `M`.
"""
function mod_prime(i, ::Val{M}) where M
    is_prime(M) || throw(DomainError(M, "modulus must be a prime number"))
    i = i % M
    return i + ifelse(signbit(i), M, 0)
end
mod_prime(i, ::Val{2}) = i & 1

"""
    Mod{M} <: Integer

Representation of finite field ``\\mathbb{Z}_M``, integers modulo small, prime `M`. Supports
integer arithmetic and can be converted to integer with `Int`.
"""
struct Mod{M} <: Integer
    value::Int

    function Mod{M}(value::Integer; check_mod=true) where M
        if check_mod
            return new{M}(mod_prime(value, Val(M)))
        else
            return new{M}(value)
        end
    end
end
Mod{M}(i::Mod{M}) where M = i

Base.Int(i::Mod) = i.value

Base.show(io::IO, i::Mod{M}) where M = print(io, Int(i), " mod ", M)

for op in (:+, :-, :*)
    @eval (Base.$op)(i::Mod{M}, j::Mod{M}) where M = Mod{M}($op(Int(i), Int(j)))
end

Base.:/(i::Mod{M}, j::Mod{M}) where M = i * inv(j)
Base.:-(i::Mod{M}) where M = Mod{M}(M - Int(i), check_mod=false)
Base.zero(::Type{Mod{M}}) where M = Mod{M}(0, check_mod=false)
Base.one(::Type{Mod{M}}) where M = Mod{M}(1, check_mod=false)
Base.sign(i::M) where M<:Mod = ifelse(iszero(i), zero(M), one(M))

Base.promote_rule(::Type{Mod{M}}, ::Type{<:Integer}) where {M} = Mod{M}

# Idea: precompute inverses and generate a function with the inverses hard-coded.
@generated function Base.inv(i::Mod{M}) where M
    err_check = quote
        iszero(i) && throw(DivideError())
    end
    if M != 2
        inverse_arr = fill(Mod{M}(0), M-1)
        inverse_arr[1] = Mod{M}(1)
        for i in 2:M-1
            inverse_arr[i] = Mod{M}(invmod(i, M))
        end
        inverse = (inverse_arr...,)

        return quote
            $err_check
            @inbounds $inverse[i]
        end
    else
        return quote
            $err_check
            i
        end
    end
end