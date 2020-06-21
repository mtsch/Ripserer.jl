"""
    IndexedSimplex{D, T, I<:Integer} <: AbstractSimplex{D, T, I}

A refinement of [`AbstractSimplex`](@ref). An indexed simplex is represented by its
dimension, diameter and combinatorial index of type `I`. It does not need to hold
information about the vertices it includes, since they can be recomputed from the index
and dimension.

By defining the [`index`](@ref), a default implementation of `sign`, `isless`,
[`vertices`](@ref) and [`coboundary`](@ref) is provided.

# Interface

* `IndexedSimplex{D[, T, I]}(index::I, diam::T)` - constructor.
* `IndexedSimplex{D[, T, I]}(vertices::NTuple{D+1, I}, diam::T)` - constructor.
* [`diam(::AbstractSimplex)`](@ref)
* [`coface_type(::AbstractSimplex)`](@ref)
* [`index(::IndexedSimplex)`](@ref)
"""
abstract type IndexedSimplex{D, T, I<:Integer} <: AbstractSimplex{D, T, I} end

"""
    index(simplex::IndexedSimplex)

Get the combinatorial index of the `simplex`. A negative index represents a simplex with
negative orientation.

```jldoctest
index(Simplex{2}((3, 2, 1), 3.2))

# output

1
```
"""
index(::IndexedSimplex)

Base.sign(sx::IndexedSimplex) = sign(index(sx))
Base.abs(sx::S) where S<:IndexedSimplex = S(abs(index(sx)), diam(sx))
Base.:-(sx::S) where S<:IndexedSimplex = S(-index(sx), diam(sx))

function Base.isless(sx1::S, sx2::S) where S<:IndexedSimplex
    return ifelse(
        diam(sx1) ≠ diam(sx2),
        diam(sx1) < diam(sx2),
        abs(index(sx1)) > abs(index(sx2))
    )
end

function Base.:(==)(sx1::IndexedSimplex{D}, sx2::IndexedSimplex{D}) where D
    return (index(sx1) == index(sx2)) & (diam(sx1) == diam(sx2))
end
function Base.isequal(sx1::IndexedSimplex{D}, sx2::IndexedSimplex{D}) where D
    return (index(sx1) == index(sx2)) & (diam(sx1) == diam(sx2))
end
function Base.hash(sx::IndexedSimplex, h::UInt64)
    return hash(index(sx), hash(diam(sx), h))
end

Base.eltype(sx::IndexedSimplex{<:Any, <:Any, I}) where I = I
Base.getindex(sx::IndexedSimplex, i) = vertices(sx)[i]
Base.firstindex(sx::IndexedSimplex) = 1
Base.lastindex(sx::IndexedSimplex{D}) where D = D + 1
Base.size(sx::IndexedSimplex{D}) where D = (D + 1,)

# vertices and indices =================================================================== #
"""
    small_binomial(n, ::Val{k})

Binomial coefficients for small, statically known values of `k`, where `n` and `k` are
always positive.
"""
small_binomial(::I, ::Val{0}) where I = one(I)
small_binomial(n, ::Val{1}) = n
function small_binomial(n::I, ::Val{k}) where {k, I}
    x = nn = I(n - k + 1)
    nn += one(I)
    for rr in I(2):I(k)
        x = div(x * nn, rr)
        nn += one(I)
    end
    return I(x)
end

"""
    find_max_vertex(idx, ::Val{k})

Use binary search to find index of first vertex in `(k-1)`-dimensional simplex with index
`idx`.
"""
function find_max_vertex(idx::I, ::Val{k}) where {I, k}
    lo = I(k - 1)
    hi = I(k + 100)
    while small_binomial(hi, Val(k)) ≤ idx
        lo = hi
        hi <<= 1
    end
    return find_max_vertex(idx, Val(k), hi + one(I), lo)
end

function find_max_vertex(idx::I, ::Val{k}, hi::I, lo::I=I(k-1)) where {I, k}
    while lo < hi - one(I)
        m = lo + ((hi - lo) >>> 0x01)
        if small_binomial(m, Val(k)) ≤ idx
            lo = m
        else
            hi = m
        end
    end
    return lo
end

"""
    vertices(index::I, ::Val{N})

Get the vertices of simplex represented by index. Returns `NTuple{N, I}`.
For regular simplices, `N` should be equal to `dim+1`!
"""
@generated function vertices(index::I, ::Val{N})::SVector{N, I} where {I, N}
    # Generate code of the form
    # index = abs(index) - 1
    # vk   = find_max_vertex(index, Val(k))
    # vk-1 = find_max_vertex(index, Val(k-1), vk)
    # ...
    # v1 = find_max_vertex(index, Val(3), v2)
    # v0 = find_max_vertex(index, Val(3), v1)
    # (vk, ..., v0) .+ 1
    vars = Symbol[Symbol("v", k) for k in N-1:-1:0]
    expr = quote
        index = abs(index) - one(I)
        $(vars[1]) = find_max_vertex(index, Val($N))
        index -= small_binomial($(vars[1]), Val($N))
    end

    for (i, k) in enumerate(N-1:-1:1)
        expr = quote
            $expr
            $(vars[i+1]) = find_max_vertex(index, Val($k), $(vars[i]))
            index -= small_binomial($(vars[i+1]), Val($k))
        end
    end
    return quote
        $expr
        SVector($(vars...)) .+ I(1)
    end
end

function vertices(sx::IndexedSimplex{D, <:Any, I}) where {D, I}
    return vertices(index(sx), Val(D+1))::SVector{D+1, I}
end

"""
    index(vertices)

Calculate the index from tuple of vertices. The index is equal to

```math
(i_d, i_{d-1}, ..., 1) \\mapsto \\sum_{k=1}^{d+1} \\binom{i_k - 1}{k},
```

where ``i_k`` are the simplex vertex indices.
"""
@generated function index(vertices::Union{NTuple{k}, SVector{k}}) where k
    # generate code of the form
    # 1 + small_binomial(vertices[1] - 1, Val(k))
    #   + small_binomial(vertices[2] - 1, Val(k-1))
    #   ...
    #   + small_binomial(vertices[k] - 1, Val(1))
    I = eltype(vertices)
    expr = quote
        one($I) + small_binomial(vertices[1] - one($I), Val($k))
    end
    for i in 2:k
        expr = quote
            $expr + small_binomial(vertices[$i] - one($I), Val($(k - i + 1)))
        end
    end
    return expr
end

index(vertex::Integer) = vertex
index(vertex::CartesianIndex) = index(vertex.I)

# (co)boundaries ========================================================================= #
struct IndexedCobounary{all_cofaces, D, I, F, S<:IndexedSimplex}
    filtration::F
    simplex::S
    vertices::NTuple{D, I}

    function IndexedCobounary{A}(
        filtration::F, simplex::S
    ) where {A, D, I, F, S<:IndexedSimplex{D, <:Any, I}}
        return new{A, D + 1, I, F, S}(filtration, simplex, Tuple(vertices(simplex)))
    end
end

function coboundary(filtration, simplex::IndexedSimplex, ::Val{A}=Val(true)) where A
    return IndexedCobounary{A}(filtration, simplex)
end

function Base.iterate(
    ci::IndexedCobounary{all_cofaces, D, I}, (v, k)=(I(n_vertices(ci.filtration) + 1), D),
) where {all_cofaces, D, I}
    @inbounds while true
        v -= one(I)
        while k > 0 && v == ci.vertices[end + 1 - k]
            all_cofaces || return nothing
            v -= one(I)
            k -= 1
        end
        v > 0 || return nothing
        sign = ifelse(iseven(k), one(I), -one(I))
        vertices = TupleTools.insertafter(ci.vertices, D - k, (v,))
        sx = simplex(ci.filtration, Val(D), vertices, sign)
        if !isnothing(sx)
            _sx::simplex_type(ci.filtration, D) = sx
            return _sx, (v, k)
        end
    end
end

struct IndexedBoundary{D, I, F, S<:IndexedSimplex}
    filtration::F
    simplex::S
    vertices::NTuple{D, I}

    function IndexedBoundary(
        filtration::F, simplex::S
    ) where {D, I, F, S<:IndexedSimplex{D, <:Any, I}}
        return new{D + 1, I, F, S}(filtration, simplex, Tuple(vertices(simplex)))
    end
end

function boundary(filtration, simplex::IndexedSimplex)
    return IndexedBoundary(filtration, simplex)
end

function Base.iterate(bi::IndexedBoundary{D, I}, k=1) where {D, I}
    while k ≤ D
        face_vertices = TupleTools.deleteat(bi.vertices, k)
        k += 1
        sign = ifelse(iseven(k), one(I), -one(I))
        sx = simplex(bi.filtration, Val(D - 2), face_vertices, sign)
        if !isnothing(sx)
            _sx::simplex_type(bi.filtration, D - 2) = sx
            return _sx, k
        end
    end
    return nothing
end

"""
    Simplex{D, T, I} <: IndexedSimplex{D, T, I}

The vanilla simplex type represented by dimension `D` and index of type `I` and a diameter
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
struct Simplex{D, T, I} <: IndexedSimplex{D, T, I}
    index ::I
    diam  ::T

    function Simplex{D, T, I}(index::Integer, diam) where {D, T, I<:Integer}
        D ≥ 0 || throw(DomainError(D, "dimension must be a non-negative integer"))
        return new{D, T, I}(I(index), T(diam))
    end
end

function Simplex{D}(index::I, diam::T) where {D, T, I<:Integer}
    return Simplex{D, T, I}(index, diam)
end
function Simplex{D}(vertices, diam::T) where {D, T}
    return Simplex{D, T, eltype(vertices)}(vertices, diam)
end
function Simplex{D, T, I}(vertices, diam) where {D, T, I<:Integer}
    length(vertices) == D + 1 ||
        throw(ArgumentError("invalid number of vertices $(length(vertices))"))
    vertices_svec = sort(SVector{D + 1}(vertices), rev=true)
    return Simplex{D, T, I}(index(vertices_svec), T(diam))
end

function Base.show(io::IO, ::MIME"text/plain", sx::Simplex{D, T, I}) where {D, T, I}
    print(io, D, "-dim Simplex", (index(sx), diam(sx)))
    print(io, ":\n  $(sign(sx) == 1 ? '+' : '-')$(vertices(sx))")
end

function Base.show(io::IO, sx::Simplex{D, M}) where {D, M}
    print(io, "Simplex{$D}($(sign(sx) == 1 ? '+' : '-')$(vertices(sx)), $(diam(sx)))")
end

# Interface implementation =============================================================== #
index(sx::Simplex) = sx.index
diam(sx::Simplex) = sx.diam
coface_type(::Type{<:Simplex{D, T, I}}) where {D, T, I} = Simplex{D+1, T, I}
face_type(::Type{<:Simplex{D, T, I}}) where {D, T, I} = Simplex{D-1, T, I}
