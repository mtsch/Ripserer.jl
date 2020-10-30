"""
    abstract type AbstractChainElement{S<:AbstractSimplex, F}

Abstract type for representing an element of a chain. Create subtypes of this if there is a
more efficient way to store a simplex and coefficient than as two fields.

An `AbstractChainElement{S, F}` behaves like a `S` when doing comparisons and like a `F`
when doing arithmetic.

Coefficient values must be in a field - they must support `+`, `-`, `*`, `/`, `inv`, `zero`,
`one` and `iszero`.

# Interface

* `AbstractChainElement{S, F}(sx::S)` - default to element with `coefficient` equal to
  `sign(sx)`.
* `AbstractChainElement{S, F}(::S, coefficient)` - convert `coefficient` to `F` if
  necessary.
* `simplex(::ChainElement)`
* `coefficient(::ChainElement)`
"""
abstract type AbstractChainElement{S<:AbstractSimplex,F<:Number} end

for op in (:+, :-)
    @eval function (Base.$op)(ce1::C, ce2::C) where {C<:AbstractChainElement}
        @boundscheck ce1 == ce2 || throw(ArgumentError("simplices don't match"))
        return C(simplex(ce1), $op(coefficient(ce1), coefficient(ce2)))
    end
end
for op in (:oneunit, :zero, :-, :+)
    @eval function (Base.$op)(elem::C) where {C<:AbstractChainElement}
        return C(simplex(elem), $op(coefficient(elem)))
    end
end
function Base.:*(elem::C, x::F) where {F<:Number,C<:AbstractChainElement{<:Any,F}}
    return C(simplex(elem), coefficient(elem) * x)
end
function Base.:*(x::F, elem::C) where {F<:Number,C<:AbstractChainElement{<:Any,F}}
    return C(simplex(elem), x * coefficient(elem))
end
function Base.:/(elem::C, x::F) where {F<:Number,C<:AbstractChainElement{<:Any,F}}
    return C(simplex(elem), coefficient(elem) / x)
end
Base.one(::Type{<:AbstractChainElement{<:Any,F}}) where {F} = one(F)
Base.one(::AbstractChainElement{<:Any,F}) where {F} = one(F)
Base.iszero(elem::AbstractChainElement) = iszero(coefficient(elem))

Base.hash(elem::AbstractChainElement, u::UInt64) = hash(simplex(elem), u)
function Base.:(==)(elem1::AbstractChainElement, elem2::AbstractChainElement)
    return simplex(elem1) == simplex(elem2)
end
function Base.isless(elem1::AbstractChainElement, elem2::AbstractChainElement)
    return isless(simplex(elem1), simplex(elem2))
end

# Make chain elements useful when they come out as representatives.
birth(elem::AbstractChainElement) = birth(simplex(elem))
index(elem::AbstractChainElement) = index(simplex(elem))
Base.sign(elem::AbstractChainElement) = sign(coefficient(elem))
vertices(elem::AbstractChainElement) = vertices(simplex(elem))

simplex_type(::Type{<:AbstractChainElement{S}}) where S = S
field_type(::Type{<:AbstractChainElement{<:Any,F}}) where F = F
simplex_type(::C) where C<:AbstractChainElement = simplex_type(C)
field_type(::C) where C<:AbstractChainElement = field_type(C)

# Make chain elements kinda similar to `Pair`s
function Base.getindex(elem::AbstractChainElement, i)
    if i == 1
        return simplex(elem)
    elseif i == 2
        return coefficient(elem)
    else
        throw(BoundsError(elem, i))
    end
end

Base.eltype(::Type{<:AbstractChainElement{S,F}}) where {S,F} = Union{S,F}
Base.iterate(elem::AbstractChainElement, i=1) = i > 2 ? nothing : (elem[i], i + 1)
Base.firstindex(::AbstractChainElement) = 1
Base.lastindex(::AbstractChainElement) = 2
Base.length(::AbstractChainElement) = 2

function Base.show(io::IO, elem::AbstractChainElement)
    return print(io, simplex(elem), " => ", coefficient(elem))
end

function Base.convert(::Type{C}, simplex::S) where {S,C<:AbstractChainElement{S}}
    return C(simplex)
end
function Base.convert(::Type{S}, elem::AbstractChainElement{S}) where {S<:AbstractSimplex}
    return simplex(elem)
end
function Base.convert(
    ::Type{C}, elem::AbstractChainElement{S,F}
) where {S,F,C<:AbstractChainElement{S,F}}
    return C(simplex(elem), coefficient(elem))
end

"""
    chain_element_type(simplex, coefficient)

Get the type of `AbstractChainElement` that is to be used for packing `simplex` and
`coefficient`. Overload this function to specialize element types for (simplex, field) type
pairs.
"""
chain_element_type(::S, ::F) where {S,F} = chain_element_type(S, F)

"""
    ChainElement{S<:AbstractSimplex, F} <: AbstractChainElement{S, F}

The default, basic subtype of [`AbstractChainElement`](@ref).
"""
struct ChainElement{S<:AbstractSimplex,F} <: AbstractChainElement{S,F}
    simplex::S
    coefficient::F

    function ChainElement{S,F}(simplex::S, coefficient=one(F)) where {S,F}
        if F <: Union{Signed, Unsigned, AbstractFloat}
            error("$F is not a field! Please try a differnet field type")
        end
        return new{S,F}(abs(simplex), sign(simplex) * F(coefficient))
    end
end

simplex(elem::ChainElement) = elem.simplex
coefficient(elem::ChainElement) = elem.coefficient
chain_element_type(::Type{S}, ::Type{F}) where {S,F} = ChainElement{S,F}

"""
    PackedElement{S<:Simplex, F<:Mod} <: AbstractChainElement{S, F}

Like `ChainElement` for `Simplex` and `Mod`. Packs the coefficient value into top `M` bits
of index. Used internally by `Chain`.
"""
struct PackedElement{S<:Simplex,F<:Mod,M,U,T} <: AbstractChainElement{S,F}
    index_coef::U
    birth::T

    function PackedElement{S,F,M,U,T}(
        simplex::S, coefficient=one(F)
    ) where {M,U,T,S<:Simplex{<:Any,T},F<:Mod{M}}
        idx = simplex.index
        coef = F(coefficient) * sign(idx)
        uidx = U(abs(idx))
        index_coef = U(Int(coef)) << (sizeof(U) * 8 - n_bits(Val(M))) | uidx

        return new{S,F,M,U,T}(index_coef, birth(simplex))
    end
end

"""
    n_bits(::Val{M})

Get numer of bits needed to represent number mod `M`.
"""
@inline n_bits(::Val{M}) where {M} = 8 * sizeof(Int) - leading_zeros(M)

function simplex(pe::PackedElement{S,<:Any,M,U}) where {S,M,U}
    mask = typemax(U) << n_bits(Val(M)) >> n_bits(Val(M))
    return S(pe.index_coef & mask, pe.birth)
end
function coefficient(elem::PackedElement{<:Any,F,M,U}) where {F,U,M}
    return F(elem.index_coef >> (sizeof(U) * 8 - n_bits(Val(M))), false)
end

function chain_element_type(
    ::Type{S}, ::Type{F}
) where {M,I,T,S<:Simplex{<:Any,T,I},F<:Mod{M}}
    if n_bits(Val(M)) ≤ 8
        return PackedElement{S,F,M,unsigned(I),T}
    else
        return ChainElement{S,F}
    end
end

function index_overflow_check(
    ::Type{S}, ::Type{F}, nv, message=""
) where {T,I,F,S<:Simplex{<:Any,T,I}}
    vertices = ntuple(i -> I(nv - i + 1), length(S))
    element = chain_element_type(S, F)(S(index(I.(vertices)), zero(T)))

    index_I = index(element)
    index_big = index(BigInt.(vertices))

    if index_I ≠ index_big
        throw(OverflowError(message))
    end
end
function index_overflow_check(::Type{<:AbstractSimplex}, F, nv, message="")
    return nothing
end
