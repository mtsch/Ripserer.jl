"""
    AbstractSimplex{D, T, I<:Integer} <: AbstractCell{I}

An abstract type for representing simplices. A simplex is an `AbstractCell` with an integer
index. A `D`-simplex has `D + 1` vertices.

If a type implements the interface below, default implementations of [`boundary`](@ref),
[`coboundary`](@ref), and [`vertices`](@ref) are provided.

# Interface

* `AbstractSimplex{0}(::I, ::T)`
* `AbstractSimplex{1}(::I, ::T)`
* [`birth(::AbstractSimplex)`](@ref)
* [`index(::AbstractSimplex)`](@ref)
* [`sign(::AbstractSimplex)`](@ref)
* `Base.:-(::AbstractSimplex)`
* [`unsafe_simplex`](@ref)
* [`unsafe_cofacet`](@ref)

"""
abstract type AbstractSimplex{D,T,I} <: AbstractCell{D,T,I} end

function Base.show(io::IO, sx::AbstractSimplex{D}) where {D}
    sgn = sign(sx) == 1 ? :+ : :-
    name = nameof(typeof(sx))
    return print(io, "$sgn$name{$D}($(vertices(sx)), $(birth(sx)))")
end
function Base.show(io::IO, ::MIME"text/plain", sx::AbstractSimplex{D}) where {D}
    sgn = sign(sx) == 1 ? :+ : :-
    name = nameof(typeof(sx))
    return print(
        io,
        "$D-dimensional $name(index=$(index(sx)), birth=$(birth(sx))):\n  $sgn",
        vertices(sx),
    )
end

###
### Vertices and indices
###
"""
    _binomial(n, ::Val{K})

Binomial coefficients for small, statically known values of `K`, where `n` and `K` are
always positive. Using `Val(K)` allows the compiler to optimize away the loop.
"""
_binomial(::I, ::Val{0}) where {I} = one(I)
_binomial(n, ::Val{1}) = n
function _binomial(n::I, ::Val{K}) where {K,I}
    x = nn = I(n - K + 1)
    nn += one(I)
    for rr in I(2):I(K)
        x = div(x * nn, rr)
        nn += one(I)
    end
    return I(x)
end

"""
    _first_vertex(index, ::Val{K})

Use binary search to find index of first vertex in `(K - 1)`-dimensional simplex with index
`index`.
"""
function _first_vertex(index::I, ::Val{K}) where {I,K}
    lo = I(K - 1)
    hi = I(K + 100)
    while _binomial(hi, Val(K)) ≤ index
        lo = hi
        hi <<= 0x01
        hi + one(I) < lo && throw(OverflowError("simplex overflowed! This is a bug"))
    end
    return _first_vertex(index, Val(K), hi + one(I), lo)
end

function _first_vertex(index::I, ::Val{K}, hi::I, lo::I=I(K - 1)) where {I,K}
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

function _first_vertex(index::I, ::Val{1}, ::I=I(0), ::I=I(0)) where {I}
    return index
end

function _first_vertex(index::I, ::Val{2}, ::I=I(0), ::I=I(0)) where {I}
    # This is https://oeis.org/A002024
    return floor(I, (√(8 * index + 1) + 1)/2)
end

#=
@inline function _vertices(index::I, ::Val{K}) where {I,K}
    index -= 1
    vk = _first_vertex(index, Val(K))
    index -= _binomial(vk, Val(K))
    return tuple(vk + 1, _vertices(index, Val(K - 1), vk)...)
end
@inline function _vertices(index::I, ::Val{K}, prev) where {I,K}
    vk = _first_vertex(index, Val(K), prev)
    index -= _binomial(vk, Val(K))
    return tuple(vk + 1, _vertices(index, Val(K - 1), vk)...)
end
@inline function _vertices(index, ::Val{1}, _)
    return (index + 1,)
end
=#
"""
    _vertices(index::I, ::Val{N})

Get the vertices of simplex represented by index. Returns `SVector{N, I}`.
For regular simplices, `N` should be equal to `dim+1`.
"""
@generated function _vertices(index::I, ::Val{N})::SVector{N,I} where {I,N}
    # Generate code of the form
    # index = abs(index) - 1
    # vk   = _first_vertex(index, Val(k))
    # vk-1 = _first_vertex(index, Val(k-1), vk)
    # ...
    # v1 = _first_vertex(index, Val(2), v2)
    # v0 = _first_vertex(index, Val(1), v1)
    # (vk, ..., v0) .+ 1
    vars = Symbol[Symbol("v", k) for k in (N - 1):-1:0]
    expr = quote
        index = abs(index) - one(I)
        $(vars[1]) = _first_vertex(index, Val($N))
        index -= _binomial($(vars[1]), Val($N))
    end

    for (i, k) in enumerate((N - 1):-1:1)
        expr = quote
            $expr
            $(vars[i + 1]) = _first_vertex(index, Val($k), $(vars[i]))
            index -= _binomial($(vars[i + 1]), Val($k))
        end
    end
    return quote
        $expr
        tuple($(vars...)) .+ I(1)
    end
end

# TODO: use lispy recursion
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
@generated function index(vertices::Union{NTuple{K},SVector{K}}) where {K}
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

#=
"""
    index_overflow_check(vertices[, message])

Check if index overflows for a particular collection of vertices. Throw an error if it does.
"""
function index_overflow_check(
    vertices,
    message="simplex $vertices overflows in $(I)! " *
            "Try using a bigger index type in your filtration.",
)
    idx = index(vertices)
    idx_big = index(BigInt.(vertices))
    if idx ≠ idx_big
        throw(OverflowError(message))
    end
end
=#

vertices(sx::AbstractSimplex) = _vertices(index(sx), Val(length(sx)))
Base.length(sx::AbstractSimplex) = length(typeof(sx))
Base.length(::Type{<:AbstractSimplex{D}}) where {D} = D + 1

###
### Coboundaries
###
function coboundary(
    filtration::AbstractFiltration, simplex::AbstractSimplex, ::Val{A}=Val(true)
) where {A}
    return Coboundary{A}(filtration, simplex)
end

struct Coboundary{A,D,I,F,S}
    filtration::F
    simplex::S
    vertices::NTuple{D,I}

    function Coboundary{A}(
        filtration::F, simplex::AbstractSimplex{D,<:Any,I}
    ) where {A,D,I,F}
        S = typeof(simplex)
        return new{A,D + 1,I,F,S}(filtration, simplex, Tuple(vertices(simplex)))
    end
end

function Base.iterate(
    ci::Coboundary{A,D,I}, (v, k)=(I(nv(ci.filtration) + 1), D)
) where {A,D,I}
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

struct SparseCoboundary{A,D,I,F,S}
    filtration::F
    simplex::S
    vertices::NTuple{D,I}
    ptrs_begin::SVector{D,I}
    ptrs_end::SVector{D,I}

    function SparseCoboundary{A}(
        filtration::F, simplex::AbstractSimplex{D,<:Any,I}
    ) where {A,D,I,F}
        verts = vertices(simplex)
        adj = adjacency_matrix(filtration)
        colptr = adj.colptr
        ptrs_begin = colptr[SVector(verts) .+ 1]
        ptrs_end = colptr[SVector(verts)]
        return new{A,D + 1,I,F,typeof(simplex)}(
            filtration, simplex, verts, ptrs_begin, ptrs_end
        )
    end
end

@propagate_inbounds @inline function _next_common(
    ptrs::SVector{D}, ptrs_end::SVector{D}, rowval
) where {D}
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
    it::SparseCoboundary{A,D,I}, (ptrs, k)=(it.ptrs_begin, D)
) where {A,D,I}
    adj = adjacency_matrix(it.filtration)
    rowval = adj.rowval
    nzval = adj.nzval
    @inbounds while true
        ptrs, v = _next_common(ptrs, it.ptrs_end, rowval)
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
            new_vertices = TupleTools.insertafter(it.vertices, D - k, (v,))
            new_edges = nzval[ptrs]
            sx = unsafe_cofacet(it.filtration, it.simplex, new_vertices, v, sign, new_edges)
            if !isnothing(sx)
                _sx::simplex_type(it.filtration, D) = sx
                return _sx, (ptrs, k)
            end
        end
    end
end

###
### Boundary
###
function boundary(filtration::AbstractFiltration, simplex::AbstractSimplex)
    return Boundary(filtration, simplex)
end

struct Boundary{D,I,F,S}
    filtration::F
    simplex::S
    vertices::NTuple{D,I}

    function Boundary(filtration::F, simplex::AbstractSimplex{D,<:Any,I}) where {D,I,F}
        S = typeof(simplex)
        return new{D + 1,I,F,S}(filtration, simplex, Tuple(vertices(simplex)))
    end
end

function Base.iterate(bi::Boundary{D,I}, k=1) where {D,I}
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
