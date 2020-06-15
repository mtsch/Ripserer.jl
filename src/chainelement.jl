"""
    abstract type AbstractChainElement{S<:AbstractSimplex, F}

Abstract type for representing an element of a chain. Create subtypes of this if there is a
more efficient way to store a simplex and coefficient than as two fields.

An `AbstractChainElement{S, F}` behaves like a `S` when doing comparisons and like a `F`
when doing arithmetic.

Coefficient values must be in a field and must support `+`, `-`, `*`, `/`, `inv`, `zero`,
`one` and `iszero`.

# Interface

* `AbstractChainElement{S, F}(sx::S)` - default to `coefficient` equal to `sign(sx)`.
* `AbstractChainElement{S, F}(::S, coefficient)` - convert `coefficient` to `F` if necessary.
* `simplex(::ChainElement)`
* `coefficient(::ChainElement)`
"""
abstract type AbstractChainElement{S<:AbstractSimplex, F<:Number} end

for op in (:+, :-)
    @eval function (Base.$op)(ce1::C, ce2::C) where C<:AbstractChainElement
        @boundscheck ce1 == ce2 || throw(ArgumentError("simplices don't match"))
        C(simplex(ce1), $op(coefficient(ce1), coefficient(ce2)))
    end
end
for op in (:oneunit, :zero, :-, :+)
    @eval (Base.$op)(ce::C) where C<:AbstractChainElement =
        C(simplex(ce), $op(coefficient(ce)))
end
function Base.:*(ce::C, x::F) where {F<:Number, C<:AbstractChainElement{<:Any, F}}
    return C(simplex(ce), coefficient(ce) * x)
end
function Base.:*(x::F, ce::C) where {F<:Number, C<:AbstractChainElement{<:Any, F}}
    return C(simplex(ce), x * coefficient(ce))
end
function Base.:/(ce::C, x::F) where {F<:Number, C<:AbstractChainElement{<:Any, F}}
    return C(simplex(ce), coefficient(ce) / x)
end
Base.one(::Type{<:AbstractChainElement{<:Any, F}}) where F = one(F)
Base.one(::AbstractChainElement{<:Any, F}) where F = one(F)
Base.iszero(ce::AbstractChainElement) = iszero(coefficient(ce))

Base.hash(ce::AbstractChainElement, u::UInt64) = hash(simplex(ce), u)
function Base.:(==)(ce1::AbstractChainElement, ce2::AbstractChainElement)
    return simplex(ce1) == simplex(ce2)
end
function Base.isless(ce1::AbstractChainElement, ce2::AbstractChainElement)
    return isless(simplex(ce1), simplex(ce2))
end

# Make chain elements useful when they come out as representatives.
diam(ce::AbstractChainElement) = diam(simplex(ce))
Base.sign(ce::AbstractChainElement) = sign(coefficient(ce))
vertices(ce::AbstractChainElement) = vertices(simplex(ce))

function Base.show(io::IO, ce::AbstractChainElement)
    print(io, simplex(ce), " => ", coefficient(ce))
end

function Base.convert(::Type{C}, simplex::S) where {S, C<:AbstractChainElement{S}}
    C(simplex)
end

"""
    chain_element_type(simplex, coefficient)

Get the type of `AbstractChainElement` that is to be used for packing `simplex` and
`coefficient`.
"""
chain_element_type(simplex::S, coefficient::F) where {S, F} = chain_element_type(S, F)

"""
    ChainElement{S<:AbstractSimplex, F} <: AbstractChainElement{S, F}

The default, basic subtype of [`AbstractChainElement`](@ref).
"""
struct ChainElement{S<:AbstractSimplex, F} <: AbstractChainElement{S, F}
    simplex::S
    coefficient::F

    function ChainElement{S, F}(simplex::S, coefficient=one(F)) where {S, F}
        return new{S, F}(abs(simplex), sign(simplex) * F(coefficient))
    end
end

# Interface specification.
simplex(ce::ChainElement) = ce.simplex
coefficient(ce::ChainElement) = ce.coefficient
chain_element_type(::Type{S}, ::Type{F}) where {S, F} = ChainElement{S, F}

# TODO: assert F is less than say 8 bits?
# TODO: check that nothing overflows - maybe in `ripserer`.
struct PackedElement{S<:IndexedSimplex, F<:Mod, M, U, T} <: AbstractChainElement{S, F}
    index_coef ::U
    diam       ::T

    function PackedElement{S, F, M, U, T}(
        simplex::S, coefficient=one(F)
    ) where {M, U, T, S<:IndexedSimplex{<:Any, T}, F<:Mod{M}}
        idx = index(simplex)
        coef = F(coefficient) * sign(idx)
        uidx = U(abs(idx))
        index_coef = Int(coef) << (sizeof(U) * 8 - n_bits(M)) | uidx

        return new{S, F, M, U, T}(index_coef, diam(simplex))
    end
end

"""
    n_bits(M)
Get numer of bits needed to represent number mod `M`.
"""
@pure n_bits(M::Int) = floor(Int, log2(M-1)) + 1

function simplex(pe::PackedElement{S, <:Any, M, U}) where {S, M, U}
    mask = typemax(U) << n_bits(M) >> n_bits(M)
    return S(pe.index_coef & mask, pe.diam)
end
function coefficient(ce::PackedElement{<:Any, F, M, U}) where {F, U, M}
    return F(ce.index_coef >> (sizeof(U) * 8 - n_bits(M)), check_mod=false)
end

function chain_element_type(
    ::Type{S}, ::Type{F}
) where {M, I, T, S<:Simplex{<:Any, T, I}, F<:Mod{M}}

    if n_bits(M) ≤ 8
        return PackedElement{S, F, M, unsigned(I), T}
    else
        return ChainElement{S, F}
    end
end
