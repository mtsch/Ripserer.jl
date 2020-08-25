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
* `Base.:-(::AbstractSimplex)`

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
julia> index(Simplex{2}((3, 2, 1), 3.2))
1
```
"""
index(::AbstractSimplex)

Base.:(==)(::AbstractSimplex, ::AbstractSimplex) = false
Base.:(==)(sx1::A, sx2::A) where A<:AbstractSimplex = index(sx1) == index(sx2)
Base.isequal(sx1::A, sx2::A) where A<:AbstractSimplex = index(sx1) == index(sx2)
Base.hash(sx::AbstractSimplex, h::UInt64) = hash(index(sx), h)

"""
    birth(simplex::AbstractSimplex)

Get the birth time of `simplex`, i.e. the time it appears in a filtration.

```jldoctest
julia> birth(Simplex{2}((3, 2, 1), 3.2))
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
always positive. Using `Val(K)` allows the compiler to optimize away the loop.
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
        hi + one(I) < lo && throw(OverflowError("simplex overflowed! This is a bug"))
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
    _vertices(index::I, ::Val{N})

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
    index(vertices)

Calculate the index from tuple or static vector of vertices. The index is equal to

```math
(i_d, i_{d-1}, ..., i_1) \\mapsto \\sum_{k=1}^{d+1} \\binom{i_k - 1}{k},
```

where ``i_k`` are the simplex vertex indices.

```jldoctest
julia> index((6,2,1))
11
```
"""
@generated function index(vertices::Union{NTuple{K}, SVector{K}}) where K
    # generate code of the form
    # 1 + _binomial(vertices[1] - 1, Val(K))
    #   + _binomial(vertices[2] - 1, Val(K-1))
    #   ...
    #   + _binomial(vertices[K] - 1, Val(1))
    I = eltype(vertices)
    expr = quote
        one($I) + _binomial(vertices[1] - one($I), Val($K))
    end
    for i in 2:K
        expr = quote
            $expr + _binomial(vertices[$i] - one($I), Val($(K - i + 1)))
        end
    end
    return expr
end

"""
    index_overflow_check(vertices[, message])

Check if index overflows for a particular collection of vertices. Throw an error if it does.
"""
function index_overflow_check(
    vertices,
    message="simplex $vertices overflows in $(I)! " *
    "Try using a bigger index type in your filtration."
)
    idx = index(vertices)
    idx_big = index(BigInt.(vertices))
    if idx ≠ idx_big
        throw(OverflowError(message))
    end
end

"""
    vertices(simplex::AbstractSimplex{dim, T, I})

Get the vertices of `simplex`. Returns `SVector{length(simplex), I}`. When `index(simplex)`
is an interger, a default implementation is provided.

# Example

```jldoctest
julia> vertices(Simplex{2}((3, 2, 1), 3.2))
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

# Examples

```jldoctest
julia> sign(Simplex{2}((3, 2, 1), 3.2))
+1

julia> sign(-Simplex{2}((3, 2, 1), 3.2))
-1
```
"""
Base.sign(::AbstractSimplex)

"""
    -(simplex::AbstractSimplex)

Reverse the simplex orientation.

# Example

```jldoctest
julia> -Simplex{2}((3, 2, 1), 3.2)
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

# Examples

```jldoctest
julia> dim(Simplex{2}((3, 2, 1), 3.2))
2

julia> dim(Cube{3, Int, 4})
3

```
"""
dim(::Type{<:AbstractSimplex{D}}) where D = D
dim(sx::AbstractSimplex) = dim(typeof(sx))

Base.abs(sx::AbstractSimplex) = sign(sx) == 1 ? sx : -sx

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
    ci::Coboundary{A, D, I}, (v, k)=(I(nv(ci.filtration) + 1), D),
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

struct SparseCoboundary{A, D, I, F, S}
    filtration::F
    simplex::S
    vertices::SVector{D, I}
    ptrs_begin::SVector{D, I}
    ptrs_end::SVector{D, I}

    function SparseCoboundary{A}(
        filtration::F, simplex::AbstractSimplex{D, <:Any, I}
    ) where {A, D, I, F}
        verts = vertices(simplex)
        adj = adjacency_matrix(filtration)
        colptr = adj.colptr
        ptrs_begin = colptr[verts .+ 1]
        ptrs_end = colptr[verts]
        return new{A, D + 1, I, F, typeof(simplex)}(
            filtration, simplex, verts, ptrs_begin, ptrs_end
        )
    end
end

@propagate_inbounds @inline function next_common(
    ptrs::SVector{D}, ptrs_end::SVector{D}, rowval
) where D
    # could also indicate when m is equal to one of the points
    ptrs = ptrs .- 1
    for i in 1:D
        ptrs[i] < ptrs_end[i] && return zero(ptrs), 0
    end
    m = rowval[ptrs[2]]
    i = 1
    while true
        ptrs_i = ptrs[i]
        row = rowval[ptrs_i]
        while row > m
            ptrs -= SVector(ntuple(isequal(i), Val(D)))
            ptrs_i -= 1
            ptrs_i < ptrs_end[i] && return zero(ptrs), zero(eltype(rowval))
            row = rowval[ptrs_i]
        end
        i = ifelse(row == m, i + 1, 1)
        i > D && return ptrs, m
        m = row
    end
end

function Base.iterate(
    it::SparseCoboundary{A, D, I}, (ptrs, k)=(it.ptrs_begin, D)
) where {A, D, I}
    adj = adjacency_matrix(it.filtration)
    rowval = adj.rowval
    nzval = adj.nzval
    @inbounds while true
        ptrs, v = next_common(ptrs, it.ptrs_end, rowval)
        if iszero(first(ptrs))
            return nothing
        elseif k > 0 && v == it.vertices[end + 1 - k]
            k -= 1
        else
            while k > 0 && v < it.vertices[end + 1 - k]
                k -= 1
            end
            !A && k ≠ D && return nothing

            sign = ifelse(iseven(k), 1, -1)
            new_vertices = insert(it.vertices, D - k + 1, v)
            new_edges = nzval[ptrs]
            sx = unsafe_cofacet(
                it.filtration, it.simplex, new_vertices, v, sign, new_edges
            )
            if !isnothing(sx)
                _sx::simplex_type(it.filtration, D) = sx
                return _sx, (ptrs, k)
            end
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
julia> Simplex{2}(2, 1)
2-dim Simplex{2}(2, 1):
  +[4, 2, 1]

julia> Simplex{10}(Int128(-10), 1.0)
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
    return Simplex{D, T, eltype(vertices)}(index(vertices_svec), T(birth))
end

index(sx::Simplex) = abs(sx.index)
birth(sx::Simplex) = sx.birth
Base.sign(sx::Simplex) = sign(sx.index)
Base.:-(sx::S) where S<:Simplex = S(-sx.index, birth(sx))
