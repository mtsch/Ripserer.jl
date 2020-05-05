# abstract simplex ======================================================================= #
"""
    AbstractSimplex{D, C, T}

An abstract type for representing simplices. A simplex is represented by its dimension,
diameter, combinatorial index and coefficient value. It does not need to hold information
about its the vertices it includes, since they can be recomputed from the index and
dimension.

`D` is the dimension, `T` is the type of distance and `C` is the coefficient type. `D` is
accessible by `dim(::AbstractSimplex)`.

# Interface

* [`index(::AbstractSimplex)`](@ref)
* [`coef(::AbstractSimplex)`](@ref)
* [`set_coef(::AbstractSimplex, ::Any)`](@ref)
* [`diam(::AbstractSimplex)`](@ref)
* [`coface_type(::AbstractSimplex)`](@ref)
* [`vertices(::AbstractSimplex)`](@ref) - optional, comes with a default implementation.
* [`coboundary(filtration, ::AbstractSimplex)`](@ref) - optional, comes with a default
  implementation.
"""
abstract type AbstractSimplex{D, C, T} end

"""
    coef(simplex::AbstractSimplex)

Get the coefficient value of `simplex`.
"""
coef(::AbstractSimplex)

"""
    set_coef(simplex::AbstractSimplex, value)

Return new `simplex` of the same type with new coefficient `value`.
"""
set_coef(::AbstractSimplex, ::Any)

"""
    diam(simplex::AbstractSimplex)

Get the diameter of `simplex`.
"""
diam(::AbstractSimplex)

"""
    index(simplex::AbstractSimplex)

Get the combinatorial index of `simplex`.
"""
index(::AbstractSimplex)

"""
    coface_type(::AbstractSimplex)
    coface_type(::Type{<:AbstractSimplex})

Get the type of coface a simplex hax. For a `D`-dimensional simplex, this is usually its
`D+1`-dimensional counterpart. Only the `coface_type(::Type{<:AbstractSimplex})` method
needs to be implemented.
"""
coface_type(sx::AbstractSimplex) =
    coface_type(typeof(sx))

"""
    dim(::AbstractSimplex)
    dim(::Type{<:AbstractSimplex})

Get the dimension of simplex i.e. the value of `D`.
"""
dim(::Type{<:AbstractSimplex{D}}) where D =
    D
dim(sx::AbstractSimplex) =
    dim(typeof(sx))

# simplex arithmetic ===================================================================== #
Base.isless(sx1::A, sx2::A) where A<:AbstractSimplex =
    diam(sx1) < diam(sx2) || diam(sx1) == diam(sx2) && index(sx1) > index(sx2)

Base.:+(sx1::A, sx2::A) where A<:AbstractSimplex =
    set_coef(sx1, coef(sx1) + coef(sx2))
Base.:-(sx1::A, sx2::A) where A<:AbstractSimplex =
    set_coef(sx1, coef(sx1) - coef(sx2))
Base.:*(sx::AbstractSimplex, x::Number) =
    set_coef(sx, coef(sx) * x)
Base.:*(x::Number, sx::AbstractSimplex) =
    set_coef(sx, x::Number * coef(sx))
Base.:-(sx::AbstractSimplex) =
    set_coef(sx, -coef(sx))
Base.:/(sx::AbstractSimplex{<:Any, C}, x::Number) where C =
    set_coef(sx, coef(sx) * inv(C(x)))
Base.one(::Type{<:AbstractSimplex{<:Any, C}}) where C =
    one(C)

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
    # vk   = find_max_vertex(index, Val(k))
    # vk-1 = find_max_vertex(index, Val(k-1), vk)
    # ...
    # v1 = find_max_vertex(index, Val(3), v2)
    # v0 = find_max_vertex(index, Val(3), v1)
    # (vk, ..., v0) .+ 1
    vars = Symbol[Symbol("v", k) for k in dim:-1:0]
    expr = quote
        index -= 1
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

"""
    vertices(sx::AbstractSimplex{dim})

Get the vertices of simplex `sx`. Returns `NTuple{dim+1, Int}`.
"""
vertices(sx::AbstractSimplex{D}) where D =
    vertices(index(sx), Val(D))::NTuple{D+1, Int}

"""
    index(vertices)

Calculate the index of vertices. The index is equal to

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

"""
    CoboundaryIterator{A, D, F, S<:AbstractSimplex{D}}

Iterator that evaluates the coboundary of a `D`-dimensional simplex. Uses the filtration to
determine which simplices are valid cofaces in the filtration. If the type parameter `A` is
`true`, return all cofaces, otherwise only return cofaces where the new vertex index is
larger than all vertex indices in `simplex`. `A=false` is used with sparse filtrations
during the call to `assemble_columns!` to generate the list of all simplices.

# Fields

* `filtration ::F`
* `simplex    ::S`
* `vertices   ::NTuple{D+1, Int}`

# Constructor

    coboundary(filtration, simplex[, Val(false)])
"""
struct CoboundaryIterator{A, D, F, S<:AbstractSimplex{D}, D2}
    filtration ::F
    simplex    ::S
    vertices   ::NTuple{D2, Int}

    CoboundaryIterator{A, D, F, S}(filtration, simplex, vertices) where {A, D, F, S} =
        new{A, D, F, S, D+1}(filtration, simplex, vertices)
end

"""
    coboundary(filtration, simplex)

Find the coboundary of `simplex`. Use the `filtration` to determine the diameters and
validity of cofaces. Iterates values of the type `coface_type(simplex)`.
"""
coboundary(
    filtration,
    simplex::S,
) where {D, M, T, S<:AbstractSimplex{D, M, T}} =
    CoboundaryIterator{true, D, typeof(filtration), S}(
        filtration, simplex, vertices(simplex)
    )
coboundary(
    filtration,
    simplex::S,
    ::Val{false},
) where {D, M, T, S<:AbstractSimplex{D, M, T}} =
    CoboundaryIterator{false, D, typeof(filtration), S}(
        filtration, simplex, vertices(simplex)
    )

function Base.iterate(ci::CoboundaryIterator{A, D},
                      (v, k)=(n_vertices(ci.filtration)+1, D+1),
                      ) where {A, D}
    diameter = ∞
    @inbounds while diameter == ∞ && v > 0
        v -= 1
        while v > 0 && v in ci.vertices
            A || return nothing
            v -= 1
            k -= 1
        end
        v == 0 && break
        diameter = diam(ci.filtration, ci.simplex, ci.vertices, v)
    end
    if diameter != ∞
        coefficient = ifelse(k % 2 == 1, -coef(ci.simplex), coef(ci.simplex))
        new_index = index(TupleTools.insertafter(ci.vertices, D+1-k, (v,)))

        coface_type(typeof(ci.simplex))(diameter, new_index, coefficient), (v, k)
    else
        nothing
    end
end
