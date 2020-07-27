"""
    AbstractSimplex{D, T, I} <: AbstractVector{I}

An abstract type for representing simplices. A simplex must have a diameter of type `T`,
which is its birth time. The dimension must be encoded in the type as `D` and can be
accessed by [`dim`](@ref).

The simplex is expected to act like an array of indices of type `I`, but this is not
actually needed for the main algorithm.

# Interface

* `AbstractSimplex{D}(::SVector{<:Any, <:Integer}, ::T)`
* [`birth(::AbstractSimplex)`](@ref)
* [`index(::AbstractSimplex)`](@ref)
* [`sign(::AbstractSimplex)`](@ref)
* [`coboundary(::Any, ::AbstractSimplex)`](@ref)
* [`boundary(::Any, ::AbstractSimplex)`](@ref)
* [`vertices(::AbstractSimplex)`](@ref)
* Base.:-(::AbstractSimplex)

"""
abstract type AbstractSimplex{D, T, I} <: AbstractVector{I} end

function Base.show(io::IO, sx::AbstractSimplex{D}) where D
    print(io, sign(sx) == 1 ? :+ : :-, nameof(typeof(sx)), "{", D, "}(",
          vertices(sx), ", ", birth(sx), ")")
end
function Base.show(io::IO, ::MIME"text/plain", sx::AbstractSimplex{D}) where D
    print(io, D, "-dimensional ", nameof(typeof(sx)),
          "(index=", index(sx), ", birth=", birth(sx), "):\n  ",
          sign(sx) == 1 ? :+ : :-, vertices(sx))
end

"""
    index(simplex::AbstractSimplex)

Get the combinatorial index of the `simplex`. The index can be any type, but should uniquely
identify a simplex. It is also used to break ties when comparing simplices with the same
birth time.

```jldoctest
index(Simplex{2}((3, 2, 1), 3.2))

# output

1
```
"""
index(::AbstractSimplex)

Base.:(==)(::AbstractSimplex, ::AbstractSimplex) = false
Base.:(==)(sx1::A, sx2::A) where A<:AbstractSimplex = index(sx1) == index(sx2)
Base.isequal(sx1::A, sx2::A) where A<:AbstractSimplex = sx1 == sx2
Base.hash(sx::AbstractSimplex, h::UInt64) = hash(index(sx), h)

"""
    birth(simplex::AbstractSimplex)

Get the birth time of `simplex`, i.e. the time it appears in a filtration.

```jldoctest
birth(Simplex{2}((3, 2, 1), 3.2))

# output

3.2
```
"""
birth(::AbstractSimplex)

function Base.isless(sx1::A, sx2::A) where {A<:AbstractSimplex}
    return ifelse(
        birth(sx1) ≠ birth(sx2),
        birth(sx1) < birth(sx2),
        index(sx1) > index(sx2),
    )
end

# vertices and indices =================================================================== #
"""
    _binomial(n, ::Val{K})

Binomial coefficients for small, statically known values of `K`, where `n` and `K` are
always positive.
"""
_binomial(::I, ::Val{0}) where I = one(I)
_binomial(n, ::Val{1}) = n
function _binomial(n::I, ::Val{K}) where {K, I}
    x = nn = I(n - K + 1)
    nn += one(I)
    for rr in I(2):I(K)
        x = div(x * nn, rr)
        nn += one(I)
    end
    return I(x)
end

"""
    _find_max_vertex(index, ::Val{K})

Use binary search to find index of first vertex in `(K - 1)`-dimensional simplex with index
`index`.
"""
function _find_max_vertex(index::I, ::Val{K}) where {I, K}
    lo = I(K - 1)
    hi = I(K + 100)
    while _binomial(hi, Val(K)) ≤ index
        lo = hi
        hi <<= 0x01
    end
    return _find_max_vertex(index, Val(K), hi + one(I), lo)
end

function _find_max_vertex(index::I, ::Val{K}, hi::I, lo::I=I(K-1)) where {I, K}
    while lo < hi - one(I)
        m = lo + ((hi - lo) >>> 0x01)
        if _binomial(m, Val(K)) ≤ index
            lo = m
        else
            hi = m
        end
    end
    return lo
end

"""
    vertices(index::I, ::Val{N})

Get the vertices of simplex represented by index. Returns `SVector{N, I}`.
For regular simplices, `N` should be equal to `dim+1`.
"""
@generated function _vertices(index::I, ::Val{N})::SVector{N, I} where {I, N}
    # Generate code of the form
    # index = abs(index) - 1
    # vk   = _find_max_vertex(index, Val(k))
    # vk-1 = _find_max_vertex(index, Val(k-1), vk)
    # ...
    # v1 = _find_max_vertex(index, Val(3), v2)
    # v0 = _find_max_vertex(index, Val(3), v1)
    # (vk, ..., v0) .+ 1
    vars = Symbol[Symbol("v", k) for k in N-1:-1:0]
    expr = quote
        index = abs(index) - one(I)
        $(vars[1]) = _find_max_vertex(index, Val($N))
        index -= _binomial($(vars[1]), Val($N))
    end

    for (i, k) in enumerate(N-1:-1:1)
        expr = quote
            $expr
            $(vars[i+1]) = _find_max_vertex(index, Val($k), $(vars[i]))
            index -= _binomial($(vars[i+1]), Val($k))
        end
    end
    return quote
        $expr
        SVector($(vars...)) .+ I(1)
    end
end

"""
    _index(vertices)

Calculate the index from tuple of vertices. The index is equal to

```math
(i_d, i_{d-1}, ..., 1) \\mapsto \\sum_{k=1}^{d+1} \\binom{i_k - 1}{k},
```

where ``i_k`` are the simplex vertex indices.
"""
@generated function _index(vertices::Union{NTuple{k}, SVector{k}}) where k
    # generate code of the form
    # 1 + _binomial(vertices[1] - 1, Val(k))
    #   + _binomial(vertices[2] - 1, Val(k-1))
    #   ...
    #   + _binomial(vertices[k] - 1, Val(1))
    I = eltype(vertices)
    expr = quote
        one($I) + _binomial(vertices[1] - one($I), Val($k))
    end
    for i in 2:k
        expr = quote
            $expr + _binomial(vertices[$i] - one($I), Val($(k - i + 1)))
        end
    end
    return expr
end

"""
    vertices(simplex::AbstractSimplex{dim, T, I})

Get the vertices of `simplex`. Returns `SVector{length(simplex), I}`. When `index(simplex)`
is an interger, a default implementation is provided.

```jldoctest
vertices(Simplex{2}((3, 2, 1), 3.2))

# output

3-element StaticArrays.SArray{Tuple{3},Int64,1,3} with indices SOneTo(3):
 3
 2
 1

```
"""
vertices(sx::AbstractSimplex) = _vertices(index(sx), Val(length(sx)))

Base.iterate(sx::AbstractSimplex) = iterate(vertices(sx))
Base.getindex(sx::AbstractSimplex, i) = vertices(sx)[i]
Base.firstindex(::AbstractSimplex) = 1
Base.lastindex(sx::AbstractSimplex) = length(sx)
Base.size(sx::AbstractSimplex) = (length(sx),)

Base.length(sx::AbstractSimplex) = length(typeof(sx))
Base.length(::Type{<:AbstractSimplex{D}}) where D = D + 1

"""
    sign(simplex::AbstractSimplex)

Get the orientation of `simplex`. Should return -1 or 1.

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
    coboundary(filtration, simplex[, Val{all_cofacets}])

Iterate over the coboundary of `simplex` in decreasing combinatorial order. Use the
`filtration` to determine the diameters and validity of cofacets.

If `all_cofacets` is `false`, only return cofaces with vertices added to the beginning of
vertex list. The method with `all_cofacets` only has to be implemented if the filtration
does not overload [`columns_to_reduce`](@ref).

Comes with a default implementation.

# Warning

If cofacets are not returned in decreasing combinatorial order, the algorithm will not work
correctly.

# Example

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

"""
    Simplex{D, T, I<:Integer} <: AbstractSimplex{D, T, I}

The vanilla simplex type represented by dimension `D` and index of type `I` and a birth time
of type `T`.

# Constructor

    Simplex{D[, T, I]}(::I, ::T)

# Examples

```jldoctest
Simplex{2}(2, 1)

# output

2-dim Simplex{2}(2, 1):
  +[4, 2, 1]
```
```jldoctest
Simplex{10}(Int128(-10), 1.0)

# output

4-dim Simplex{3}(1.0, 10, 2):
  -Int128[12, 11, 10, 9, 8, 7, 6, 5, 4, 2, 1]
```
"""
struct Simplex{D, T, I<:Integer} <: AbstractSimplex{D, T, I}
    index::I
    birth::T

    function Simplex{D, T, I}(index::Integer, birth) where {D, T, I<:Integer}
        D ≥ 0 || throw(DomainError(D, "dimension must be a non-negative integer"))
        return new{D, T, I}(I(index), convert(T, birth))
    end
end

function Simplex{D}(index::I, birth::T) where {D, T, I<:Integer}
    return Simplex{D, T, I}(index, birth)
end
function Simplex{D}(vertices, birth::T) where {D, T}
    length(vertices) == D + 1 ||
        throw(ArgumentError("invalid number of vertices $(length(vertices))"))
    vertices_svec = sort(SVector{D + 1}(vertices), rev=true)
    return Simplex{D, T, eltype(vertices)}(_index(vertices_svec), T(birth))
end

index(sx::Simplex) = abs(sx.index)
birth(sx::Simplex) = sx.birth
Base.sign(sx::Simplex) = sign(sx.index)
Base.:-(sx::S) where S<:Simplex = S(-sx.index, birth(sx))
