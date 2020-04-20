# distance matrix stuff ================================================================== #
"""
    edge_lt(e1, e2) =

Compare edges like DiameterSimplexComparer.

* by increasing diameter,
* by decreasing combinatorial index.
"""
edge_lt((e1, i1), (e2, i2)) =
    e1 < e2 || e1 == e2 && i1 > i2

function edges(dist::AbstractMatrix{T}, thresh=typemax(T)) where T
    n = size(dist, 1)
    res = Tuple{T, Tuple{Int, Int}}[]
    @inbounds for j in 1:n, i in j+1:n
        l = dist[i, j]
        l ≤ thresh && push!(res, (l, (i, j)))
    end
    sort!(res, lt=edge_lt)
end

function edges(dist::AbstractSparseMatrix{T}, thresh=typemax(T)) where T
    res = Tuple{T, Tuple{Int, Int}}[]
    I, J, V = findnz(dist)
    for (i, j, l) in zip(I, J, V)
        i > j || continue
        l ≤ thresh && push!(res, (l, (i, j)))
    end
    sort!(res, lt=edge_lt)
end

"""
    is_distance_matrix(dist)

Return true if dist is a valid distance matrix.
"""
is_distance_matrix(dist) =
    issymmetric(dist) && all(iszero(dist[i, i]) for i in 1:size(dist, 1))

"""
    distances(metric, points)

Return distance matrix calculated from `points` with `metric`.
"""
function distances(metric, points)
    dim = length(first(points))
    T = eltype(first(points))
    pairwise(metric, reshape(reinterpret(T, points), (dim, length(points))), dims=2)
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
