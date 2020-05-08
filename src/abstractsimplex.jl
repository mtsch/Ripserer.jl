"""
    AbstractSimplex{D, T}

An abstract type for representing simplices. A simplex must have a diameter of type `T`,
which is its birth time. The dimension must be encoded in the type as `D` and can be
accessed by [`dim`](@ref).

# Interface

* `AbstractSimplex{D}(::NTuple{D+1, <:Integer}, ::T)`
* [`diam(::AbstractSimplex)`](@ref)
* [`Base.sign(::AbstractSimplex)`](@ref)
* [`Base.:-(::AbstractSimplex)`](@ref)
* [`coface_type(::AbstractSimplex)`](@ref)
* [`vertices(::AbstractSimplex)`](@ref)
* [`coboundary(::Any, ::AbstractSimplex)`](@ref)
* `Base.isless(::AbstractSimplex, ::AbstractSimplex)`
"""
abstract type AbstractSimplex{D, T} end

"""
    diam(simplex::AbstractSimplex)

Get the diameter of `simplex`.
"""
diam(::AbstractSimplex)

"""
    sign(simplex::AbstractSimplex)

Get the orientation of `simplex`. Returns -1 or 1.
"""
Base.sign(::AbstractSimplex)

"""
    -(simplex::AbstractSimplex)

Reverse the simplex orientation.
"""
Base.:-(::AbstractSimplex)
Base.:+(sx::AbstractSimplex) = sx

Base.:(==)(::AbstractSimplex, ::AbstractSimplex) = false
Base.:(==)(sx1::AbstractSimplex{D}, sx2::AbstractSimplex{D}) where D =
    diam(sx1) == diam(sx2) && vertices(sx1) == vertices(sx2)
Base.hash(sx::AbstractSimplex, h::UInt64) =
    hash(vertices(sx), hash(diam(sx), h))

"""
    coface_type(::AbstractSimplex)
    coface_type(::Type{<:AbstractSimplex})

Get the type of coface a simplex hax. For a `D`-dimensional simplex, this is usually its
`D+1`-dimensional counterpart. Only the method for the type needs to be implemented.
"""
coface_type(sx::AbstractSimplex) = coface_type(typeof(sx))

"""
    dim(::AbstractSimplex)
    dim(::Type{<:AbstractSimplex})

Get the dimension of simplex i.e. the value of `D`.
"""
dim(::Type{<:AbstractSimplex{D}}) where D = D
dim(sx::AbstractSimplex) = dim(typeof(sx))

Base.abs(sx::AbstractSimplex) = sign(sx) == 1 ? sx : -sx
