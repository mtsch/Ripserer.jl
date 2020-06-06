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
* `Base.isless(::AbstractSimplex, ::AbstractSimplex)`
* [`coface_type(::AbstractSimplex)`](@ref)
* [`vertices(::AbstractSimplex)`](@ref)
* [`coboundary(::Any, ::AbstractSimplex)`](@ref)
"""
abstract type AbstractSimplex{D, T} end

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
  -(3, 2, 1)
```
"""
Base.:-(::AbstractSimplex)
Base.:+(sx::AbstractSimplex) = sx

Base.:(==)(::AbstractSimplex, ::AbstractSimplex) = false
function Base.:(==)(sx1::AbstractSimplex{D}, sx2::AbstractSimplex{D}) where D
    return diam(sx1) == diam(sx2) && vertices(sx1) == vertices(sx2)
end
Base.hash(sx::AbstractSimplex, h::UInt64) = hash(vertices(sx), hash(diam(sx), h))

"""
    coface_type(::AbstractSimplex)
    coface_type(::Type{<:AbstractSimplex})

Get the type of simplex's coface. For a `D`-dimensional simplex, this is usually its
`D+1`-dimensional counterpart. Only the method for the type needs to be implemented.

```jldoctest
coface_type(Simplex{2}((3, 2, 1), 3.2))

# output

Simplex{3, Float64, Int}
```
"""
coface_type(sx::AbstractSimplex) = coface_type(typeof(sx))

"""
    face_type(::AbstractSimplex)
    face_type(::Type{<:AbstractSimplex})

Get the type of a simplex's face. For a `D`-dimensional simplex, this is usually its
`D-1`-dimensional counterpart. Only the method for the type needs to be implemented.

```jldoctest
coface_type(Simplex{2}((3, 2, 1), 3.2))

# output

Simplex{1, Float64, Int}
```
"""
face_type(sx::AbstractSimplex) = face_type(typeof(sx))

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

Get the vertices of `simplex`. Returns `NTuple{dim+1, Int}`. In the algorithm, only the
method for 2-simplices is actually used.

```jldoctest
vertices(Simplex{2}((3, 2, 1), 3.2))

# output

(3, 2, 1)
```
"""
vertices(::AbstractSimplex)


"""
    coboundary(filtration, simplex[, Val{all_cofaces}])

Iterate over the coboundary of `simplex`. Use the `filtration` to determine the diameters
and validity of cofaces. Iterates values of the type `coface_type(simplex)`. If
`all_cofaces` is `false`, only return cofaces with vertices added to the beginning of vertex
list.

```jldoctest coboundary
filtration = RipsFiltration([0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0])

for c in coboundary(filtration, Simplex{1}(2, 1))
    println(c)
end

# output

Simplex{2}(+(4, 3, 1), 1)
Simplex{2}(-(3, 2, 1), 1)
```

```jldoctest coboundary
for c in coboundary(filtration, Simplex{1}(2, 1), Val(false))
    println(c)
end

# output

Simplex{2}(+(4, 3, 1), 1)
```
"""
coboundary(::Any, ::AbstractSimplex)
