"""
    AbstractSimplex{D, T, I} <: AbstractVector{I}

An abstract type for representing simplices. A simplex must have a diameter of type `T`,
which is its birth time. The dimension must be encoded in the type as `D` and can be
accessed by [`dim`](@ref).

The simplex is expected to act like an array of indices of type `I`, but this is not
actually needed for the main algorithm.

# Interface

* `AbstractSimplex{D}(::SVector{<:Any, <:Integer}, ::T)`
* [`diam(::AbstractSimplex)`](@ref)
* [`Base.sign(::AbstractSimplex)`](@ref)
* [`Base.:-(::AbstractSimplex)`](@ref)
* `Base.isless(::AbstractSimplex, ::AbstractSimplex)`
* [`vertices(::AbstractSimplex)`](@ref)
* [`coboundary(::Any, ::AbstractSimplex)`](@ref)
* [`boundary(::Any, ::AbstractSimplex)`](@ref)
* `length(::Type{AbstractSimplex})` - only if `D`-simplex does not have `D + 1` vertices.
"""
abstract type AbstractSimplex{D, T, I} <: AbstractVector{I} end

Base.:(==)(::AbstractSimplex, ::AbstractSimplex) = false
function Base.:(==)(sx1::A, sx2::A) where A<:AbstractSimplex
    return diam(sx1) == diam(sx2) && vertices(sx1) == vertices(sx2)
end
function Base.isequal(sx1::A, sx2::A) where A<:AbstractSimplex
    return diam(sx1) == diam(sx2) && vertices(sx1) == vertices(sx2)
end
Base.hash(sx::AbstractSimplex, h::UInt64) = hash(vertices(sx), hash(diam(sx), h))

Base.iterate(sx::AbstractSimplex) = iterate(vertices(sx))
Base.getindex(sx::AbstractSimplex, i) = vertices(sx)[i]
Base.firstindex(::AbstractSimplex) = 1
Base.lastindex(sx::AbstractSimplex) = length(sx)
Base.size(sx::AbstractSimplex) = (length(sx),)

Base.length(sx::AbstractSimplex) = length(typeof(sx))
Base.length(::Type{<:AbstractSimplex{D}}) where D = D + 1

function Base.show(io::IO, sx::AbstractSimplex{D}) where D
    print(io, nameof(typeof(sx)), "{", D, "}(",
          sign(sx) == 1 ? :+ : :-, vertices(sx), ", ", diam(sx), ")")
end

"""
    diam(simplex::AbstractSimplex)

Get the diameter of `simplex`.

```jldoctest
diam(Simplex{2}((3, 2, 1), 3.2))

# output

3.2
```
"""
diam(::AbstractSimplex)

"""
    sign(simplex::AbstractSimplex)

Get the orientation of `simplex`. Returns -1 or 1.

```jldoctest
sign(Simplex{2}((3, 2, 1), 3.2))

# output

+1
```
"""
Base.sign(::AbstractSimplex)

"""
    -(simplex::AbstractSimplex)

Reverse the simplex orientation.

```jldoctest
-Simplex{2}((3, 2, 1), 3.2)

# output

2-dim Simplex(1, 1):
  -[3, 2, 1]
```
"""
Base.:-(::AbstractSimplex)
Base.:+(sx::AbstractSimplex) = sx

"""
    dim(::AbstractSimplex)
    dim(::Type{<:AbstractSimplex})

Get the dimension of simplex i.e. the value of `D`.

```jldoctest
dim(Simplex{2}((3, 2, 1), 3.2))

# output

2
```
"""
dim(::Type{<:AbstractSimplex{D}}) where D = D
dim(sx::AbstractSimplex) = dim(typeof(sx))

Base.abs(sx::AbstractSimplex) = sign(sx) == 1 ? sx : -sx

"""
    vertices(simplex::AbstractSimplex{dim})

Get the vertices of `simplex`. Returns `SVector{length(simplex), Int}`.

```jldoctest
vertices(Simplex{2}((3, 2, 1), 3.2))

# output

3-element StaticArrays.SArray{Tuple{3},Int64,1,3} with indices SOneTo(3):
 3
 2
 1

```
"""
vertices(::AbstractSimplex)

"""
    coboundary(filtration, simplex[, Val{all_cofacets}])

Iterate over the coboundary of `simplex` in decreasing combinatorial order. Use the
`filtration` to determine the diameters and validity of cofacets. If `all_cofacets` is
`false`, only return cofaces with vertices added to the beginning of vertex list.

Comes with a default implementation.

```jldoctest coboundary
filtration = Rips([0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0])

for c in coboundary(filtration, Simplex{1}(2, 1))
    println(c)
end

# output

Simplex{2}(+[4, 3, 1], 1)
Simplex{2}(-[3, 2, 1], 1)
```

```jldoctest coboundary
for c in coboundary(filtration, Simplex{1}(2, 1), Val(false))
    println(c)
end

# output

2-dim Simplex(1, 1):
  +[4, 3, 1]
```
"""
function coboundary(filtration, simplex::AbstractSimplex, ::Val{A}=Val(true)) where A
    return Coboundary{A}(filtration, simplex)
end

"""
    boundary(filtration, simplex[, Val{all_cofacets}])

Iterate over the boundary of `simplex`. Use the `filtration` to determine the diameters
and validity of cofacets.

Comes with a default implementation.

```jldoctest boundary
filtration = Rips([0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0])

for f in boundary(filtration, Simplex{2}(2, 1))
    println(f)
end

# output

Simplex{2}(+[2, 1], 1)
Simplex{2}(-[4, 1], 1)
Simplex{2}(+[4, 2], 1)
```
"""
boundary(filtration, simplex::AbstractSimplex) = Boundary(filtration, simplex)

# (co)boundaries ========================================================================= #
struct Coboundary{A, D, I, F, S}
    filtration::F
    simplex::S
    vertices::NTuple{D, I}

    function Coboundary{A}(
        filtration::F, simplex::AbstractSimplex{D, <:Any, I}
    ) where {A, D, I, F}
        S = typeof(simplex)
        return new{A, D + 1, I, F, S}(filtration, simplex, Tuple(vertices(simplex)))
    end
end

function Base.iterate(
    ci::Coboundary{A, D, I}, (v, k)=(I(n_vertices(ci.filtration) + 1), D),
) where {A, D, I}
    @inbounds while true
        v -= one(I)
        while k > 0 && v == ci.vertices[end + 1 - k]
            A || return nothing
            v -= one(I)
            k -= 1
        end
        v > 0 || return nothing
        sign = ifelse(iseven(k), one(I), -one(I))
        new_vertices = TupleTools.insertafter(ci.vertices, D - k, (v,))
        sx = unsafe_cofacet(ci.filtration, ci.simplex, new_vertices, v, sign)
        if !isnothing(sx)
            _sx::simplex_type(ci.filtration, D) = sx
            return _sx, (v, k)
        end
    end
end

struct Boundary{D, I, F, S}
    filtration::F
    simplex::S
    vertices::NTuple{D, I}

    function Boundary(
        filtration::F, simplex::AbstractSimplex{D, <:Any, I}
    ) where {D, I, F}
        S = typeof(simplex)
        return new{D + 1, I, F, S}(filtration, simplex, Tuple(vertices(simplex)))
    end
end

function Base.iterate(bi::Boundary{D, I}, k=1) where {D, I}
    while k ≤ D
        facet_vertices = TupleTools.deleteat(bi.vertices, k)
        k += 1
        sign = ifelse(iseven(k), one(I), -one(I))
        sx = unsafe_simplex(bi.filtration, Val(D - 2), facet_vertices, sign)
        if !isnothing(sx)
            _sx::simplex_type(bi.filtration, D - 2) = sx
            return _sx, k
        end
    end
    return nothing
end
