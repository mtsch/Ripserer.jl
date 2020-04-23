# Beware: black magic ahead.
using Ripserer
using Ripserer: Coboundary, dist_type
import Ripserer: index, vertices, diam, coef, set_coef

using TupleTools

# prime fields =========================================================================== #
"""
    is_prime(n)

Return `true` if `n` is a prime number.
"""
function is_prime(n)
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

"""
    prime_check(::Val{M})

Check that type parameter `M` is prime at compile time.
"""
@generated function prime_check(::Val{M}) where M
    is_prime(M) || throw(DomainError(M, "modulus not prime"))
    nothing
end

"""
    mod_prime(i, ::Val{M})

Like `mod`, but with prime `M`.
"""
function mod_prime(i, ::Val{M}) where M
    prime_check(Val(M))
    i = i % M
    i + ifelse(signbit(i), M, 0)
end
mod_prime(i, ::Val{2}) =
    i & 1

"""
    PrimeField{M} <: Integer

Representation of finite field, integers modulo `M`.
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
Base.promote_rule(::Type{PrimeField{M}}, ::Type{Int}) where {D, M} =
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

# magic simplex ========================================================================== #
"""
    n_bits(M)

Get numer of bits needed to represent number mod `M`.
"""
n_bits(M) =
    floor(Int, log2(M-1)) + 1

"""
    MagicSimplex{D, M, T} <: AbstractSimplex{PrimeField{M}, T}

The vanilla simplex type with coefficient values from ``\\mathbb{Z}_M``, integers modulo
prime `M`. `index` and `coef` are packed into a single `UInt64`.

# Constructor

    MagicSimplex{D, M}(::T, index::Integer, coef::Integer)
"""
struct MagicSimplex{D, M, T} <: AbstractSimplex{PrimeField{M}, T}
    diam       ::T
    index_coef ::UInt64

    @generated function MagicSimplex{D, M, T}(
        diam::T, index::Integer, coef::Integer
    ) where {D, M, T}

        prime_check(Val(M))
        bits = n_bits(M)
        Expr(:new, :(MagicSimplex{D, M, T}),
             :diam,
             :(UInt64(index) << $bits + mod_prime(coef, Val(M))))
    end
    @generated function MagicSimplex{D, M, T}(
        diam::T, index::Integer, coef::PrimeField{M}
    ) where {D, M, T}

        bits = n_bits(M)
        Expr(:new, :(MagicSimplex{D, M, T}),
             :diam,
             :(UInt64(index) << $bits + Int(coef)))
    end
end

MagicSimplex{D, M}(diam::T, index, coef) where {D, M, T} =
    MagicSimplex{D, M, T}(diam, index, coef)

@generated function index(sx::MagicSimplex{<:Any, M}) where M
    shift = n_bits(M)
    :(reinterpret(Int64, sx.index_coef >> $shift))
end

@generated function coef(sx::MagicSimplex{<:Any, M}) where M
    mask = 1 << n_bits(M) - 1
    :(PrimeField{$M}(sx.index_coef & $mask, check_mod=false))
end

diam(sx::MagicSimplex) =
    sx.diam

set_coef(sx::MagicSimplex{D, M, T}, coef) where {D, M, T} =
    MagicSimplex{D, M, T}(diam(sx), index(sx), coef)

Base.show(io::IO, sx::MagicSimplex{D, M}) where {D, M} =
    print(io, "MagicSimplex{", D, ", ", M, "}", (diam(sx), index(sx), Int(coef(sx))))

# magical stuff ========================================================================== #
magic_binomial(_, ::Val{0}) = 1
magic_binomial(n, ::Val{1}) = n

@generated function magic_binomial(n, ::Val{k}) where k
    coefficients = [0]
    k_minus_1 = k - 1

    quote
        n_coef = length($coefficients)
        if n_coef < n + 1
            resize!($coefficients, n)
            for i in n_coef+1:n
                $coefficients[i] = magic_binomial(i-1, Val($k_minus_1)) + $coefficients[i-1]
            end
        end
        @inbounds $coefficients[n]
    end
end

# Generate code of the form
# vk   = find_max_vertex(index, Val(k))
# vk-1 = find_max_vertex(index, Val(k-1), vk)
# ...
# v1 = find_max_vertex(index, Val(3), v2)
# v0 = find_max_vertex(index, Val(3), v1)
# (vk, ..., v0)
@generated function vertices(index, ::Val{dim}) where dim
    vars = Symbol[Symbol("v", k) for k in dim:-1:0]
    expr = quote
        index -= 1
        $(vars[1]) = find_max_vertex(index, Val($dim+1))
        index -= magic_binomial($(vars[1]), Val($dim+1))
    end

    for (i, k) in enumerate(dim:-1:1)
        expr = quote
            $expr
            $(vars[i+1]) = find_max_vertex(index, Val($k), $(vars[i]) - 1)
            index -= magic_binomial($(vars[i+1]), Val($k))
        end
    end
    quote
        $expr
        tuple($(vars...)) .+ 1
    end
end

vertices(sx::MagicSimplex{D}) where D =
    vertices(index(sx), Val(D))::NTuple{D+1, Int}

function find_max_vertex(idx, ::Val{k}) where k
    n_max = 10
    while magic_binomial(n_max, Val(k)) ≤ idx
        n_max <<= 1
    end
    find_max_vertex(idx, Val(k), n_max)
end

function find_max_vertex(idx, ::Val{k}, n_max) where k
    hi = n_max + 1
    lo = k - 1
    while lo < hi - 1
        m = lo + ((hi - lo) >>> 0x01)
        if magic_binomial(m, Val(k)) ≤ idx
            lo = m
        else
            hi = m
        end
    end
    lo
end

# generate code of the form
# 1 + magic_binomial(vertices[1] - 1, Val(k))
#   + magic_binomial(vertices[2] - 1, Val(k-1))
#   ...
#   + magic_binomial(vertices[k] - 1, Val(1))
@generated function index(vertices::NTuple{k}) where k
    expr = quote
        1 + magic_binomial(vertices[1]-1, Val($k))
    end
    for i in 2:k
        expr = quote
            $expr + magic_binomial(vertices[$i]-1, Val($(k - i + 1)))
        end
    end
    expr
end

struct MagicCoboundary{D1, D2, M, T, F, S<:MagicSimplex{D1, M, T}}
    filtration ::F
    simplex    ::S
    vertices   ::NTuple{D2, Int}
end

coboundary(filtration, simplex::MagicSimplex{dim, M}) where {dim, M} =
    MagicCoboundary{dim, dim+1, M, dist_type(filtration), typeof(filtration), typeof(simplex)}(
        filtration,
        simplex,
        vertices(simplex)
    )

function Base.iterate(ci::MagicCoboundary{D1, D2, M, T}, (v, k)=(n_vertices(ci.filtration), D1+1)) where {D1, D2, M, T}
    diameter = ∞
    @inbounds while diameter == ∞ && v > 0
        while v > 0 && v in ci.vertices
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

        MagicSimplex{D2, M, T}(diameter, new_index, coefficient), (v - 1, k)
    else
        nothing
    end
end

Base.@pure coface_type(::Type{MagicSimplex{D, M, T}}) where {D, M, T} =
    MagicSimplex{D+1, M, T}

struct MagicCoboundary2{D, F<:AbstractFiltration, S<:MagicSimplex{D}, D2}
    filtration ::F
    simplex    ::S
    vertices   ::NTuple{D2, Int}

    function MagicCoboundary2{D, F, S}(filtration, simplex, vertices) where {D, F, S}
        new{D, F, S, D+1}(filtration, simplex, vertices)
    end
end

coboundary2(
    filtration::AbstractFiltration{T, S},
    simplex::S,
) where {D, M, T, S<:MagicSimplex{D, M, T}} =
    MagicCoboundary2{D, typeof(filtration), S}(filtration, simplex, vertices(simplex))

function Base.iterate(ci::MagicCoboundary2{D},
                      (v, k)=(n_vertices(ci.filtration), D+1),
                      ) where D
    diameter = ∞
    @inbounds while diameter == ∞ && v > 0
        while v > 0 && v in ci.vertices
            v -= 1
            k -= 1
        end
        v == 0 && break
        diameter = diam(ci.filtration, ci.simplex, ci.vertices, v)
        v -= 1
    end
    if diameter != ∞
        coefficient = ifelse(k % 2 == 1, -coef(ci.simplex), coef(ci.simplex))
        # todo: be smarter than sort...
        new_index = index(TupleTools.sort(tuple(ci.vertices..., v-1), rev=true))

        coface_type(typeof(ci.simplex))(diameter, new_index, coefficient), (v, k)
    else
        nothing
    end
end

# ======================================================================================== #
# benching
include(joinpath(@__DIR__, "../test/data.jl"))
Random.seed!(7350)

function count_cofaces(coboundary, sx)
    count = 0
    for i in 1:10000
        for coface in coboundary(sx, 2)
            count += 1
        end
    end
    count ÷ 10000
end
# Distances are between 0 and 2.
dists = rand_dist_matrix(4000)
sx = Simplex{2}(dists[1, 2], 1, 1)
msx = MagicSimplex{2, 2}(dists[1, 2], 1, 1)

coboundary_full_nothreshold = Coboundary(RipsFiltration(dists, threshold=10), 2)
coboundary_full_threshold1 = Coboundary(RipsFiltration(dists, threshold=1), 2)
coboundary_sparse_75 = Coboundary(SparseRipsFiltration(dists, threshold=1.5), 2)
coboundary_sparse_50 = Coboundary(SparseRipsFiltration(dists, threshold=1), 2)
coboundary_sparse_25 = Coboundary(SparseRipsFiltration(dists, threshold=0.5), 2)

filtration_full_nothreshold = RipsFiltration(dists, threshold=10, simplex_type=typeof(msx))
filtration_full_threshold1 = RipsFiltration(dists, threshold=1, simplex_type=typeof(msx))
filtration_sparse_75 = SparseRipsFiltration(dists, threshold=1.5, simplex_type=typeof(msx))
filtration_sparse_50 = SparseRipsFiltration(dists, threshold=1, simplex_type=typeof(msx))
filtration_sparse_25 = SparseRipsFiltration(dists, threshold=0.5, simplex_type=typeof(msx))

function show_cofaces()
    sx = MagicSimplex{2, 2}(3.0, 1, 1)
    flt = RipsFiltration(icosahedron, simplex_type=typeof(sx))
    for coface in coboundary(flt, sx)
        @show vertices(coface)
    end
end

function count_cofaces_new(filtration, sx)
    count = 0
    for i in 1:10000
        for coface in coboundary(filtration, sx)
            count += 1
        end
    end
    count ÷ 10000
end

function count_cofaces_newer(filtration, sx)
    count = 0
    for i in 1:10000
        for coface in coboundary2(filtration, sx)
            count += 1
        end
    end
    count ÷ 10000
end

suite_old = BenchmarkGroup()
suite_new = BenchmarkGroup()

suite_old["full, no threshold"] =
    @benchmarkable count_cofaces($coboundary_full_nothreshold, $sx)
suite_old["full, threshold=1"] =
    @benchmarkable count_cofaces($coboundary_full_threshold1, $sx)
suite_old["sparse, 75% full"] =
    @benchmarkable count_cofaces($coboundary_sparse_75, $sx)
suite_old["sparse, 50% full"] =
    @benchmarkable count_cofaces($coboundary_sparse_50, $sx)
suite_old["sparse, 25% full"] =
    @benchmarkable count_cofaces($coboundary_sparse_25, $sx)

suite_new["full, no threshold"] =
    @benchmarkable count_cofaces_newer($filtration_full_nothreshold, $msx)
suite_new["full, threshold=1"] =
    @benchmarkable count_cofaces_newer($filtration_full_threshold1, $msx)
suite_new["sparse, 75% full"] =
    @benchmarkable count_cofaces_newer($filtration_sparse_75, $msx)
suite_new["sparse, 50% full"] =
    @benchmarkable count_cofaces_newer($filtration_sparse_50, $msx)
suite_new["sparse, 25% full"] =
    @benchmarkable count_cofaces_newer($filtration_sparse_25, $msx)

#=
println("tuning new")
println(tune!(suite_new))
println("tuning old")
println(tune!(suite_old))
=#

println("run old")
b_old = run(suite_old)
println(b_old)
println("run new")
b_new = run(suite_new)
println(b_new)
