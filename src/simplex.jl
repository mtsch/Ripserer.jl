# prime fields =========================================================================== #
"""
    is_prime(n)

Return `true` if `n` is a prime number.
"""
@pure function is_prime(n::Int)
    if iseven(n) || n < 2
        n == 2
    else
        p = 3
        q = n / p
        while p ≤ q
            iszero(n % p) && return false
            p += 2
            q = n / p
        end
        true
    end
end
@pure is_prime(::Any) =
    false

"""
    mod_prime(i, ::Val{M})

Like `mod`, but with prime `M`.
"""
function mod_prime(i, ::Val{M}) where M
    is_prime(M) || throw(DomainError(M, "modulus must be a prime number"))
    i = i % M
    i + ifelse(signbit(i), M, 0)
end
mod_prime(i, ::Val{2}) =
    i & 1

"""
    PrimeField{M} <: Integer

Representation of finite field ``\\mathbb{Z}_M``, integers modulo `M`. Supports integer
arithmetic and can be converted to integer with `Int`.
"""
struct PrimeField{M} <: Integer
    value::Int

    function PrimeField{M}(value::Integer; check_mod=true) where M
        if check_mod
            new{M}(mod_prime(value, Val(M)))
        else
            new{M}(value)
        end
    end
end
PrimeField{M}(i::PrimeField{M}) where M =
    i

Base.Int(i::PrimeField) =
    i.value

Base.show(io::IO, i::PrimeField{M}) where M =
    print(io, Int(i), " mod ", M)

for op in (:+, :-, :*)
    @eval (Base.$op)(i::PrimeField{M}, j::PrimeField{M}) where M =
        PrimeField{M}($op(Int(i), Int(j)))
end

Base.:/(i::PrimeField{M}, j::PrimeField{M}) where M =
    i * inv(j)
Base.:-(i::PrimeField{M}) where M =
    PrimeField{M}(M - Int(i), check_mod=false)
Base.zero(::Type{PrimeField{M}}) where M =
    PrimeField{M}(0, check_mod=false)
Base.one(::Type{PrimeField{M}}) where M =
    PrimeField{M}(1, check_mod=false)
Base.promote_rule(::Type{PrimeField{M}}, ::Type{Int}) where {M} =
    PrimeField{M}

# Idea: precompute inverses and generate a function with the inverses hard-coded.
@generated function Base.inv(i::PrimeField{M}) where M
    err_check = quote
        iszero(i) && throw(DivideError())
    end
    if M != 2
        inverse_arr = fill(PrimeField{M}(0), M-1)
        inverse_arr[1] = PrimeField{M}(1)
        for i in 2:M-1
            inverse_arr[i] = PrimeField{M}(invmod(i, M))
        end
        inverse = (inverse_arr...,)

        quote
            $err_check
            @inbounds $inverse[i]
        end
    else
        quote
            $err_check
            i
        end
    end
end

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
    vertices(sx::AbstractSimplex{dim})

Get the vertices of simplex represented by index. Returns `NTuple{dim, Int}`.
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

Get the vertices of simplex `sx`. Returns `NTuple{dim, Int}`.
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
    filtration::AbstractFiltration{T},
    simplex::S,
) where {D, M, T, S<:AbstractSimplex{D, M, T}} =
    CoboundaryIterator{true, D, typeof(filtration), S}(
        filtration, simplex, vertices(simplex)
    )
coboundary(
    filtration::AbstractFiltration{T},
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
        # todo: be smarter than sort...
        new_index = index(TupleTools.sort(tuple(ci.vertices..., v), rev=true))

        coface_type(typeof(ci.simplex))(diameter, new_index, coefficient), (v, k)
    else
        nothing
    end
end


# default simplex ======================================================================== #
"""
    n_bits(M)

Get numer of bits needed to represent number mod `M`.
"""
@pure n_bits(M::Int) =
    floor(Int, log2(M-1)) + 1

"""
    Simplex{D, M, T} <: AbstractSimplex{D, PrimeField{M}, T}

The vanilla simplex type with coefficient values from ``\\mathbb{Z}_M``, integers modulo
`M`, where `M` is prime. `index` and `coef` are packed into a single `U<:Unsigned`.

# Constructor

    Simplex{D, M}(::T, index::Integer, coef::Integer)

# Examples

```jldoctest
julia> Simplex{2, 2}(1, 2, 3)
2-dim Simplex{2}(1, 2, 1):
  (4, 2, 1)

julia> Simplex{10, 3}(1.0, Int128(10), 2)
4-dim Simplex{3}(1.0, 10, 2) with UInt128 index:
  (12, 11, 10, 9, 8, 7, 6, 5, 4, 2, 1)
`˙`
"""
struct Simplex{D, M, T, U<:Unsigned} <: AbstractSimplex{D, PrimeField{M}, T}
    diam       ::T
    index_coef ::U

    function Simplex{D, M, T, U}(
        diam::T, index::I, coef::Integer
    ) where {D, M, T, I<:Integer, U}
        is_prime(M) || throw(DomainError(M, "modulus must be a prime number"))
        D > 0 || throw(DomainError(D, "dimension must be a positive integer"))
        U === unsigned(I) || throw(DomainError(U, "type parameters must match"))
        bits = n_bits(M)
        new{D, M, T, unsigned(I)}(diam, unsigned(index) << bits + mod_prime(coef, Val(M)))
    end

    function Simplex{D, M, T, U}(
        diam::T, index::I, coef::PrimeField{M}
    ) where {D, M, T, I<:Integer, U}
        D > 0 || throw(DomainError(D, "dimension must be a positive integer"))
        U === unsigned(I) || throw(DomainError(U, "type parameters must match"))
        bits = n_bits(M)
        new{D, M, T, unsigned(I)}(diam, unsigned(index) << bits + I(coef))
    end
end

Simplex{D, M}(diam::T, index::Integer, coef::Integer) where {D, M, T} =
    Simplex{D, M, T}(diam, index, coef)
Simplex{D, M, T}(diam::T, index::I, coef::Integer) where {D, M, T, I<:Integer} =
    Simplex{D, M, T, unsigned(I)}(diam, index, coef)

function Base.show(io::IO, ::MIME"text/plain", sx::Simplex{D, M, <:Any, U}) where {D, M, U}
    print(io, D, "-dim Simplex{", M, "}", (diam(sx), index(sx), Int(coef(sx))))
    if !(U <: UInt)
        print(io, " with ", U, " index")
    end
    print(io, ":\n  ", vertices(sx))
end

Base.show(io::IO, sx::Simplex{D, M}) where {D, M} =
    print(io, "Simplex{", D, ", ", M, "}",
          (diam(sx), index(sx), Int(coef(sx))), " = ", vertices(sx))

index(sx::Simplex{<:Any, M}) where M =
    signed(sx.index_coef >> n_bits(M))

function coef(sx::Simplex{<:Any, M}) where M
    mask = 1 << n_bits(M) - 1
    PrimeField{M}(sx.index_coef & mask, check_mod=false)
end

diam(sx::Simplex) =
    sx.diam

set_coef(sx::Simplex{D, M, T}, coef) where {D, M, T} =
    Simplex{D, M, T}(diam(sx), index(sx), coef)

@pure coface_type(::Type{Simplex{D, M, T, U}}) where {D, M, T, U} =
    Simplex{D+1, M, T, U}
