"""
    IndexedSimplex{D, T, I<:Integer} <: AbstractSimplex{D, T}

A refinement of [`AbstractSimplex`](@ref). An indexed simplex is represented by its
dimension, diameter and combinatorial index. It does not need to hold information about its
the vertices it includes, since they can be recomputed from the index and dimension.

By defining the [`index`](@ref), a default implementation of `sign`, `isless`,
[`vertices`](@ref) and [`coboundary`](@ref) is provided.

# Interface

* `IndexedSimplex{D[, T, I]}(index::I, diam::T)` - constructor.
* `IndexedSimplex{D[, T, I]}(vertices::NTuple{D+1, I}, diam::T)` - constructor.
* [`diam(::AbstractSimplex)`](@ref)
* [`coface_type(::AbstractSimplex)`](@ref)
* [`index(::IndexedSimplex)`](@ref)
"""
abstract type IndexedSimplex{D, T, I<:Integer} <: AbstractSimplex{D, T} end

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

Base.sign(sx::IndexedSimplex) =
    sign(index(sx))

Base.:-(sx::S) where S<:IndexedSimplex =
    S(-index(sx), diam(sx))

Base.isless(sx1::S, sx2::S) where S<:IndexedSimplex =
    ifelse(diam(sx1) ≠ diam(sx2),
           diam(sx1) < diam(sx2),
           abs(index(sx1)) > abs(index(sx2)))

Base.:(==)(sx1::IndexedSimplex{D}, sx2::IndexedSimplex{D}) where D =
    (index(sx1) == index(sx2)) & (diam(sx1) == diam(sx2))
Base.hash(sx::IndexedSimplex, h::UInt64) =
    hash(index(sx), hash(diam(sx), h))

# vertices and indices =================================================================== #
"""
    small_binomial(n, ::Val{k})

Binomial coefficients for small, statically known values of `k`, where `n` and `k` are
always positive.
"""
small_binomial(_, ::Val{0}) = 1
small_binomial(n, ::Val{1}) = n
function small_binomial(n, ::Val{k}) where k
    n0, k0 = n, k
    sgn = 1
    x = nn = n - k + 1
    nn += 1
    for rr in 2:k
        x = div(x * nn, rr)
        nn += 1
    end
    x
end

"""
    find_max_vertex(idx, ::Val{k})

Use binary search to find index of first vertex in `(k-1)`-dimensional simplex with index
`idx`.
"""
function find_max_vertex(idx, ::Val{k}) where k
    lo = k - 1
    hi = k + 100
    while small_binomial(hi, Val(k)) ≤ idx
        lo = hi
        hi <<= 1
    end
    find_max_vertex(idx, Val(k), hi + 1, lo)
end

function find_max_vertex(idx, ::Val{k}, hi, lo=k-1) where k
    while lo < hi - 1
        m = lo + ((hi - lo) >>> 0x01)
        if small_binomial(m, Val(k)) ≤ idx
            lo = m
        else
            hi = m
        end
    end
    lo
end

"""
    vertices(index, ::Val{dim})

Get the vertices of simplex represented by index. Returns `NTuple{dim+1, Int}`.
"""
@generated function vertices(index, ::Val{dim}) where dim
    # Generate code of the form
    # index = abs(index) - 1
    # vk   = find_max_vertex(index, Val(k))
    # vk-1 = find_max_vertex(index, Val(k-1), vk)
    # ...
    # v1 = find_max_vertex(index, Val(3), v2)
    # v0 = find_max_vertex(index, Val(3), v1)
    # (vk, ..., v0) .+ 1
    vars = Symbol[Symbol("v", k) for k in dim:-1:0]
    expr = quote
        index = abs(index) - 1
        $(vars[1]) = find_max_vertex(index, Val($dim+1))
        index -= small_binomial($(vars[1]), Val($dim+1))
    end

    for (i, k) in enumerate(dim:-1:1)
        expr = quote
            $expr
            $(vars[i+1]) = find_max_vertex(index, Val($k), $(vars[i]))
            index -= small_binomial($(vars[i+1]), Val($k))
        end
    end
    quote
        $expr
        tuple($(vars...)) .+ 1
    end
end

vertices(sx::IndexedSimplex{D}) where D =
    vertices(index(sx), Val(D))::NTuple{D+1, Int}

"""
    index(vertices)

Calculate the index from tuple of vertices. The index is equal to

```math
(i_d, i_{d-1}, ..., 1) \\mapsto \\sum_{k=1}^{d+1} \\binom{i_k - 1}{k},
```

where ``i_k`` are the simplex vertex indices.
"""
@generated function index(vertices::NTuple{k}) where k
    # generate code of the form
    # 1 + small_binomial(vertices[1] - 1, Val(k))
    #   + small_binomial(vertices[2] - 1, Val(k-1))
    #   ...
    #   + small_binomial(vertices[k] - 1, Val(1))
    expr = quote
        1 + small_binomial(vertices[1]-1, Val($k))
    end
    for i in 2:k
        expr = quote
            $expr + small_binomial(vertices[$i]-1, Val($(k - i + 1)))
        end
    end
    expr
end

# coboundaries =========================================================================== #
"""
    IndexedCobounary{all_cofaces, D, F, S<:IndexedSimplex{D-1}}

Iterator that evaluates the coboundary of a `(D - 1)`-dimensional indexed simplex. Uses the
filtration to determine which simplices are valid cofaces in the filtration and to determine
their diameter. If the type parameter `all_cofaces` is `true`, return all cofaces, otherwise
only return cofaces where the new vertex index is larger than all vertex indices in
`simplex`. `all_cofaces=false` is used with during the call to [`assemble_columns!`](@ref)
to generate the list of all simplices.

# Fields

* `filtration ::F`
* `simplex    ::S`
* `vertices   ::NTuple{D, Int}`

# Constructor

    coboundary(filtration, simplex[, Val(false)])
"""
struct IndexedCobounary{all_cofaces, D, F, S<:IndexedSimplex}
    filtration ::F
    simplex    ::S
    vertices   ::NTuple{D, Int}

    IndexedCobounary{A}(filtration::F, simplex::S) where {A, D, F, S<:IndexedSimplex{D}} =
        new{A, D + 1, F, S}(filtration, simplex, vertices(simplex))
end

coboundary(filtration, simplex::IndexedSimplex) =
    IndexedCobounary{true}(filtration, simplex)
coboundary(filtration, simplex::IndexedSimplex, ::Val{false}) =
    IndexedCobounary{false}(filtration, simplex)

function Base.iterate(
    ci::IndexedCobounary{all_cofaces, D}, (v, k)=(n_vertices(ci.filtration) + 1, D),
) where {all_cofaces, D}
    diameter = ∞
    @inbounds while diameter == ∞ && v > 0
        v -= 1
        while v > 0 && v in ci.vertices
            all_cofaces || return nothing
            v -= 1
            k -= 1
        end
        v == 0 && break
        diameter = diam(ci.filtration, ci.simplex, ci.vertices, v)
    end
    if diameter != ∞
        sign = ifelse(iseven(k), 1, -1)
        new_index = index(TupleTools.insertafter(ci.vertices, D - k, (v,))) * sign

        return coface_type(ci.simplex)(new_index, diameter), (v, k)
    else
        return nothing
    end
end

"""
    Simplex{D, T, I} <: IndexedSimplex{D, T, I}

The vanilla simplex type represented by dimension `D` and index of type `I` and a diameter
of type `T`.

# Constructor

    Simplex{D[, T, I]}(::I, ::T)

# Examples

```jldoctest
julia> Simplex{2}(2, 1)
2-dim Simplex{2}(2, 1):
  +(4, 2, 1)

julia> Simplex{10}(Int128(-10), 1.0)
4-dim Simplex{3}(1.0, 10, 2) with UInt128 index:
  -(12, 11, 10, 9, 8, 7, 6, 5, 4, 2, 1)
`˙`
"""
struct Simplex{D, T, I} <: IndexedSimplex{D, T, I}
    index ::I
    diam  ::T

    function Simplex{D, T, I}(index, diam) where {D, T, I}
        D ≥ 0 || throw(DomainError(D, "dimension must be a non-negative integer"))
        new{D, T, I}(I(index), T(diam))
    end
end

Simplex{D}(index::I, diam::T) where {D, T, I<:Integer} =
    Simplex{D, T, I}(index, diam)
Simplex{D}(vertices::NTuple{<:Any, I}, diam::T) where {D, T, I<:Integer} =
    Simplex{D, T, I}(vertices, diam)
function Simplex{D, T, I}(vertices::NTuple{N}, diam) where {D, T, I<:Integer, N}
    N == D + 1 || throw(ArgumentError("invalid number of vertices"))

    Simplex{D, T, I}(index(vertices), T(diam))
end

function Base.show(io::IO, ::MIME"text/plain", sx::Simplex{D, T, I}) where {D, T, I}
    print(io, D, "-dim Simplex", (index(sx), diam(sx)))
    if I ≢ Int64
        print(io, " with ", I, " index")
    end
    print(io, ":\n  $(sign(sx) == 1 ? '+' : '-')$(vertices(sx))")
end

Base.show(io::IO, sx::Simplex{D, M}) where {D, M} =
    print(io, "Simplex{$D}($(sign(sx) == 1 ? '+' : '-')$(vertices(sx)), $(diam(sx)))")

# Interface implementation =============================================================== #
index(sx::Simplex) = sx.index

diam(sx::Simplex) = sx.diam

coface_type(::Type{<:Simplex{D, T, I}}) where {D, T, I} = Simplex{D+1, T, I}
