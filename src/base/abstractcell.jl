"""
    AbstractCell{D, T, I}

An abstract type for representing simplices and other types of cells, such as cubes. A cell
is determined by the following.

* Its dimension encoded in the type parameter `D`. Can be accessed with [`dim`](@ref)

* The `index` of type `I` is used to differentiate between different cells. In general, two
  cells of the same type (and hence dimension) and the same `index` are considered to be
  equal.

* The `birth` of type `T` determines when the cell enters the filtration. Note that two
  cells with the same `index` should also have the same `birth`.

* `vertices` should return the cell's vertices as a tuple. For 0-cells (vertices), the
  `vertices` are also used to index into a filtration's
  [`vertices`](vertices(::AbstractFiltration)).

* The `sign` determens its orientation. Note that `cell == -cell`.

# Interface

* [`birth(::AbstractCell)`](@ref)::T
* [`index(::AbstractCell)`](@ref)::I
* [`vertices(σ::AbstractCell)`](@ref)::NTuple{length(σ),I}
* length(::Type{AbstractCell})
* [`sign(::AbstractCell)`](@ref)
* [`Base.:-(::AbstractCell)`](@ref)
* [`coboundary(::Any, ::AbstractCell)`](@ref)
* [`boundary(::Any, ::AbstractCell)`](@ref)

"""
abstract type AbstractCell{D,T,I} end

"""
    dim(::AbstractCell)
    dim(::Type{<:AbstractCell})

Get the dimension of a cell i.e. the value of `D`. Can also be called on the type.

# Examples

```jldoctest
julia> dim(Simplex{2}((3, 2, 1), 3.2))
2

julia> dim(Cube{3, Int, 4})
3

```
"""
dim(σ::AbstractCell) = dim(typeof(σ))
dim(::Type{<:AbstractCell{D}}) where {D} = D

index_type(σ) = index_type(typeof(σ))
index_type(::Type{<:AbstractCell{<:Any,<:Any,I}}) where {I} = I
birth_type(σ) = birth_type(typeof(σ))
birth_type(::Type{<:AbstractCell{<:Any,T}}) where {T} = T

"""
    index(σ::AbstractCell)

Get the combinatorial index of the `σ`. The index can be any type, but should uniquely
identify a cell. It is also used to break ties when comparing simplices with the same birth
time.

```jldoctest
julia> index(Simplex{2}((3, 2, 1), 3.2))
1
```
"""
index(::AbstractCell)

"""
    birth(σ::AbstractCell)

Get the birth time of `σ`, i.e. the time it first appears in the filtration.

# Example

```jldoctest
julia> birth(Simplex{2}((3, 2, 1), 3.2))
3.2
```
"""
birth(::AbstractCell)

"""
    vertices(σ::AbstractCell{D,T,I})::NTuple{length(σ),I}

Get the vertices of `σ`.

# Example

```jldoctest
julia> vertices(Simplex{2}((3, 2, 1), 3.2))
(3, 2, 1)

```
"""
vertices(::AbstractCell)

Base.length(σ::AbstractCell) = length(typeof(σ))

"""
    sign(σ::AbstractCell)

Get the orientation of `σ`. Should return -1 or 1.

# Examples

```jldoctest
julia> sign(Simplex{2}((3, 2, 1), 3.2))
1

julia> sign(-Simplex{2}((3, 2, 1), 3.2))
-1
```
"""
Base.sign(::AbstractCell)

"""
    -(σ::AbstractCell)

Reverse the cell orientation.

# Example

```jldoctest
julia> -Simplex{2}((3, 2, 1), 3.2)
2-dimensional Simplex(index=1, birth=3.2):
  -(3, 2, 1)

```
"""
Base.:-(::AbstractCell)
Base.:+(σ::AbstractCell) = σ
Base.abs(σ::AbstractCell) = sign(σ) == 1 ? σ : -σ

# Equality and order
Base.:(==)(::AbstractCell, ::AbstractCell) = false
Base.:(==)(σ::A, τ::A) where {A<:AbstractCell} = index(σ) == index(τ)
Base.isequal(σ::A, τ::A) where {A<:AbstractCell} = index(σ) == index(τ)
Base.hash(σ::AbstractCell, h::UInt64) = hash(index(σ), h)
function Base.isless(σ::A, τ::A) where {A<:AbstractCell}
    return ifelse(birth(σ) ≠ birth(τ), birth(σ) < birth(τ), index(σ) > index(τ))
end

# Iteration and indexing
Base.to_index(σ::AbstractCell) = SVector(vertices(σ))
Base.eltype(::Type{S}) where {S<:AbstractCell} = index_type(S)
function Base.iterate(σ::AbstractCell, (vs, i)=(vertices(σ), 1))
    if i > length(vs)
        return nothing
    else
        return vs[i], (vs, i + 1)
    end
end

"""
    coboundary(filtration, simplex[, Val{all_cofacets}])

Iterate over the coboundary of `simplex` by decreasing `index`. Use the
`filtration` to determine the diameters and validity of cofacets.

If `all_cofacets` is `false`, only return cofaces with vertices added to the beginning of
vertex list. The method with `all_cofacets` only has to be implemented if the filtration
does not overload [`columns_to_reduce`](@ref).

Comes with a default implementation.

!!! warning
    If cofacets are not returned in decreasing `index` order, the algorithm will not work
    correctly. If there is no avoiding it, define `emergent_pairs(...) = false` for your
    filtration.

# Examples

```jldoctest coboundary
filtration = Rips([0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0])

for c in Ripserer.coboundary(filtration, Simplex{1}(2, 1))
    println(c)
end

# output

+Simplex{2}((4, 3, 1), 1)
-Simplex{2}((3, 2, 1), 1)
```

```jldoctest coboundary
for c in Ripserer.coboundary(filtration, Simplex{1}(2, 1), Val(false))
    println(c)
end

# output

+Simplex{2}((4, 3, 1), 1)
```
"""
coboundary

"""
    boundary(filtration, simplex[, Val{all_cofacets}])

Iterate over the boundary of `simplex` by increasing `index`. Use the `filtration` to
determine the diameters and validity of cofacets.

Comes with a default implementation.

!!! warning
    If facets are not returned in increasing `index` order, the (homology) algorithm will
    not work correctly. If there is no avoiding it, define `emergent_pairs(...) = false` for
    your filtration.

# Example

```jldoctest boundary
filtration = Rips([0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0])

for f in Ripserer.boundary(filtration, Simplex{2}(2, 1))
    println(f)
end

# output

+Simplex{1}((2, 1), 1)
-Simplex{1}((4, 1), 1)
+Simplex{1}((4, 2), 1)
```
"""
boundary
