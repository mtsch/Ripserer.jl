"""
    abstract type AbstractChainElement{S<:AbstractSimplex, F}

Abstract type for representing an element of a chain. Create subtypes of this if there is a
more efficient way to store a simplex and coefficient than as two fields.

An `AbstractChainElement{S, F}` behaves like a `S` when doing comparisons and like a `F`
when doing arithmetic.

Coefficient values must be in a field and must support `+`, `-`, `*`, `/`, `inv`, `zero`,
`one` and `iszero`.

# Interface

* `AbstractChainElement{S, F}(sx::S)` - default to `coef` equal to `sign(sx)`.
* `AbstractChainElement{S, F}(::S, coef)` - convert `coef` to `F` if necessary.
* `simplex(::ChainElement)`
* `coef(::ChainElement)`
"""
abstract type AbstractChainElement{S<:AbstractSimplex, F<:Number} end

for op in (:+, :-)
    @eval (Base.$op)(ce1::C, ce2::C) where C<:AbstractChainElement =
        C(simplex(ce1), $op(coef(ce1), coef(ce2)))
end
for op in (:oneunit, :zero, :-, :+)
    @eval (Base.$op)(ce::C) where C<:AbstractChainElement =
        C(simplex(ce), $op(coef(ce)))
end
Base.:*(ce::C, x::F) where {F<:Number, C<:AbstractChainElement{<:Any, F}} =
    C(simplex(ce), coef(ce) * x)
Base.:*(x::F, ce::C) where {F<:Number, C<:AbstractChainElement{<:Any, F}} =
    C(simplex(ce), x * coef(ce))
Base.:/(ce::C, x::F) where {F<:Number, C<:AbstractChainElement{<:Any, F}} =
    C(simplex(ce), coef(ce) / x)
Base.one(::Type{<:AbstractChainElement{<:Any, F}}) where F =
    one(F)
Base.one(::AbstractChainElement{<:Any, F}) where F =
    one(F)
Base.iszero(ce::AbstractChainElement) =
    iszero(coef(ce))

Base.hash(ce::AbstractChainElement, u::UInt64) =
    hash(simplex(ce), u)
Base.:(==)(ce1::AbstractChainElement, ce2::AbstractChainElement) =
    simplex(ce1) == simplex(ce2)
Base.isless(ce1::AbstractChainElement, ce2::AbstractChainElement) =
    isless(simplex(ce1), simplex(ce2))

"""
    chain_element_type(simplex, coef)

Get the type of `AbstractChainElement` that is to be used for packing `simplex` and
`coefficient`.
"""
chain_element_type(simplex::S, coef::F) where {S, F} = chain_element_type(S, F)

"""
    ChainElement{S<:AbstractSimplex, F} <: AbstractChainElement{S, F}

The default, basic subtype of [`AbstractChainElement`](@ref).
"""
struct ChainElement{S<:AbstractSimplex, F} <: AbstractChainElement{S, F}
    simplex::S
    coef::F

    ChainElement{S, F}(simplex::S, coef) where {S, F} =
        new{S, F}(abs(simplex), sign(simplex) * F(coef))

    ChainElement{S, F}(simplex::S) where {S, F} =
        new{S, F}(abs(simplex), F(sign(simplex)))
end

simplex(ce::ChainElement) = ce.simplex
coef(ce::ChainElement) = ce.coef

chain_element_type(::Type{S}, ::Type{F}) where {S, F} =
    ChainElement{S, F}

# TODO? efficient packing of Simplex and PrimeField?
