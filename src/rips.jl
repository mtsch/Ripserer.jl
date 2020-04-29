# distance matrix stuff ================================================================== #
"""
    edges(dist::AbstractMatrix{T}, thresh, S)

Return sorted edges of type `S` in distance matrix with length lower than `thresh`.
"""
function edges(dist::AbstractMatrix{T}, thresh, S) where T
    n = size(dist, 1)
    res = S[]
    @inbounds for j in 1:n, i in j+1:n
        l = dist[i, j]
        l ≤ thresh && push!(res, S(l, index((i, j)), 1))
    end
    sort!(res)
end

function edges(dist::AbstractSparseMatrix{T}, thresh, S) where T
    res = S[]
    I, J, V = findnz(dist)
    for (i, j, l) in zip(I, J, V)
        i > j || continue
        l ≤ thresh && push!(res, S(l, index((i, j)), 1))
    end
    sort!(res)
end

"""
    distances(metric, points[, births])

Return distance matrix calculated from `points` with `metric`. If given, add birth times
from `births` to the diagonal.
"""
function distances(metric, points, births=nothing)
    dim = length(first(points))
    T = eltype(first(points))
    dists = pairwise(metric, reshape(reinterpret(T, points), (dim, length(points))), dims=2)
    if !isnothing(births)
        for i in 1:length(points)
            dists[i, i] = births[i]
        end
    end
    dists
end

# infinity =============================================================================== #
"""
    Infinity

`Infinity()` is bigger than _anything_ else, except `missing` and `Inf`. It is used to:

* Avoiding using `typemax(T)` in persistence intervals. Getting death times of
  `9223372036854775807` doesn't look good.
* Returned by `diam(::AbstractFiltration, args...)` to signal that a simplex should be
  skipped.
"""
struct Infinity end

Base.show(io::IO, ::Infinity) =
    print(io, "∞")

(::Type{T})(::Infinity) where T<:AbstractFloat =
    typemax(T)
for op in (:<, :>, :isless, :isequal, :(==))
    @eval (Base.$op)(x::Real, ::Infinity) =
        $op(x, Inf)
    @eval (Base.$op)(::Infinity, x::Real) =
        $op(Inf, x)
end
Base.isapprox(::Infinity, x::Real; args...) =
    isapprox(Inf, x; args...)
Base.isapprox(x::Real, ::Infinity; args...) =
    isapprox(x, Inf; args...)

Base.isless(::Infinity, ::Missing) =
    false
Base.isless(::Missing, ::Infinity) =
    true
Base.isless(::Infinity, ::Infinity) =
    false
Base.:>(::Infinity, ::Infinity) =
    false
Base.isless(a, ::Infinity) =
    true
Base.isless(::Infinity, a) =
    true

Base.isfinite(::Infinity) =
    false

const ∞ = Infinity()

# flag =================================================================================== #
"""
    AbstractFlagFiltration{T, S} <: AbstractFiltration{T, S}

An abstract flag filtration is a filtration of flag complexes. Its subtypes can overload
`dist(::AbstractFlagFiltration{T}, u, v)::Union{T, Infinity}` instead of `diam`.
`diam(::AbstractFlagFiltration, ...)` defaults to maximum `dist` among vertices.
"""
abstract type AbstractFlagFiltration{T, S} <: AbstractFiltration{T, S} end

@propagate_inbounds function diam(flt::AbstractFlagFiltration, vertices)
    n = length(vertices)
    res = typemin(dist_type(flt))
    for i in 1:n, j in i+1:n
        d = dist(flt, vertices[j], vertices[i])
        res = ifelse(res > d, res, d)
    end
    ifelse(res > threshold(flt), ∞, res)
end

@propagate_inbounds function diam(flt::AbstractFlagFiltration, sx::AbstractSimplex, us, v)
    res = diam(sx)
    for u in us
        # Even though this looks like a tight loop, v changes way more often than us, so
        # this is the faster order of indexing by u and v.
        d = dist(flt, v, u)
        res = ifelse(res > d, res, d)
    end
    ifelse(res > threshold(flt), ∞, res)
end

edges(flt::AbstractFlagFiltration) =
    edges(flt.dist, threshold(flt), edge_type(flt))

"""
    dist(::AbstractFlagFiltration, u, v)

Return the distance between vertices `u` and `v`. If the distance is higher than the
threshold, return `Infinity()` instead.
"""
dist(::AbstractFlagFiltration, ::Any, ::Any)

# rips =================================================================================== #
"""
    default_rips_threshold(dists)

The default threshold is equal to the radius of the input space. At this threshold, all
vertices are connected to a vertex `x` and the homology becomes trivial. If any of the
distances is negative, default threshold defaults to `typemax(eltype(dists))`.
"""
default_rips_threshold(dists::AbstractMatrix{T}) where T =
    minimum(maximum(abs, dists[:, i]) for i in 1:size(dists, 1))

"""
    RipsFiltration{T, S<:AbstractSimplex{<:Any, T}} <: AbstractFlagFiltration{T, S}

This type represents a filtration of Vietoris-Rips complexes.
Diagonal items are treated as vertex birth times.

# Constructor

    RipsFiltration(
        distance_matrix;
        modulus=2,
        threshold=default_rips_threshold(dist),
        edge_type=Simplex{1, modulus, T, Int, UInt}
    )
"""
struct RipsFiltration{
    T, S<:AbstractSimplex{1, <:Any, T}, A<:AbstractMatrix{T}
} <: AbstractFlagFiltration{T, S}

    dist      ::A
    threshold ::T
end

function RipsFiltration(
    dist::AbstractMatrix{T};
    modulus=2,
    threshold=default_rips_threshold(dist),
    edge_type::DataType=Simplex{1, modulus, T, UInt64}
) where T

    issymmetric(dist) ||
        throw(ArgumentError("`dist` must be a distance matrix"))
    is_prime(modulus) ||
        throw(ArgumentError("`modulus` must be prime"))
    edge_type <: AbstractSimplex{1, <:Any, T} ||
        throw(ArgumentError("`edge_type` must be a subtype of `AbstractSimplex{1}`"))
    !issparse(dist) ||
        throw(ArgumentError("`dits` is sparse. Use `SparseRipsFiltration` instead"))
    RipsFiltration{T, edge_type, typeof(dist)}(dist, T(threshold))
end

n_vertices(rips::RipsFiltration) =
    size(rips.dist, 1)

@propagate_inbounds dist(rips::RipsFiltration{T}, i, j) where T =
    ifelse(i == j, zero(T), rips.dist[i, j])

threshold(rips::RipsFiltration) =
    rips.threshold

birth(rips::RipsFiltration, i) =
    rips.dist[i, i]

"""
    SparseRipsFiltration{T, S<:AbstractSimplex{<:Any, T}} <: AbstractFlagFiltration{T, S}

This type represents a filtration of Vietoris-Rips complexes.
The distance matrix will be converted to a sparse matrix with all values greater than
threshold deleted. Off-diagonal zeros in the matrix are treated as ∞. Diagonal items are
treated as vertex birth times.

# Constructor

    SparseRipsFiltration(
        distance_matrix;
        modulus=2,
        threshold=nothing,
        eltype=Simplex{modulus, T}
    )
"""
struct SparseRipsFiltration{
    T, S<:AbstractSimplex{1, <:Any, T}, A<:AbstractSparseMatrix{T}
}<: AbstractFlagFiltration{T, S}

    dist ::A
end

function SparseRipsFiltration(
    dist::AbstractMatrix{T};
    modulus=2,
    threshold=nothing,
    edge_type::DataType=Simplex{1, modulus, T, UInt64},
) where T
    issymmetric(dist) ||
        throw(ArgumentError("`dist` must be a distance matrix"))
    is_prime(modulus) ||
        throw(ArgumentError("`modulus` must be prime"))
    edge_type <: AbstractSimplex{1, <:Any, T} ||
        throw(ArgumentError("`edge_type` must be a subtype of `AbstractSimplex{1}`"))

    new_dist = sparse(dist)
    if !isnothing(threshold)
        births = diag(dist)
        SparseArrays.fkeep!(new_dist, (_, _ , v) -> v ≤ threshold)
        for i in 1:size(dist, 1)
            new_dist[i, i] = births[i]
        end
    end
    SparseRipsFiltration{T, edge_type, typeof(new_dist)}(new_dist)
end

n_vertices(rips::SparseRipsFiltration) =
    size(rips.dist, 1)

@propagate_inbounds function dist(rips::SparseRipsFiltration{T}, i, j) where T
    res = rips.dist[i, j]
    ifelse(i == j, zero(T), ifelse(iszero(res), ∞, res))
end

# Threshold was handled by deleting entries in the matrix.
threshold(rips::SparseRipsFiltration) =
    ∞

SparseArrays.issparse(::Type{<:SparseRipsFiltration}) =
    true

@propagate_inbounds function diam(rips::SparseRipsFiltration, sx::AbstractSimplex, us, v)
    res = diam(sx)
    for u in us
        # Since indexing in sparse matrices is expensive, we want to abort the loop early
        # even though the number of vertices in us is small.
        d = dist(rips, v, u)
        d == ∞ && return ∞
        res = ifelse(res > d, res, d)
    end
    res
end

birth(rips::SparseRipsFiltration, i) =
    rips.dist[i, i]
