"""
    AbstractRipsFiltration{I<:Signed, T} <: AbstractFiltration{I, T}

An abstract Vietoris-Rips filtration. Its subtypes can only overload [`dist`](@ref) and get
default implementations for the rest of the filtration interface.
"""
abstract type AbstractRipsFiltration{I<:Signed, T} <: AbstractFiltration{I, T} end

"""
    dist(::AbstractRipsFiltration, u, v)
    dist(::AbstractRipsFiltration)

Return the distance between vertices `u` and `v`. If the distance is somehow invalid, it may
return `missing` instead. If `u` and `v` are not given, return the distance matrix.
"""
dist(::AbstractRipsFiltration, ::Any, ::Any)

n_vertices(rips::AbstractRipsFiltration) = size(dist(rips), 1)

birth(rips::AbstractRipsFiltration) = diag(dist(rips))

simplex_type(::Type{<:AbstractRipsFiltration{I, T}}, D) where {I, T} = Simplex{D, T, I}

@inline @propagate_inbounds function unsafe_simplex(
    ::Type{S}, rips::AbstractRipsFiltration{I, T}, vertices, sign
) where {I, T, S<:Simplex}
    if dim(S) == 0
        return S(I(sign) * vertices[1], birth(rips, vertices[1]))
    else
        n = length(vertices)
        diameter = typemin(T)
        for i in 1:n, j in i+1:n
            d = dist(rips, vertices[j], vertices[i])
            if ismissing(d) || d > threshold(rips)
                return nothing
            else
                _d::T = d
                diameter = ifelse(_d > diameter, _d, diameter)
            end
        end
        return S(I(sign) * index(vertices), diameter)
    end
end

@inline @propagate_inbounds function unsafe_cofacet(
    ::Type{S},
    rips::AbstractRipsFiltration{I, T},
    simplex,
    cofacet_vertices,
    new_vertex,
    sign,
) where {I, T, S<:Simplex}
    diameter = birth(simplex)
    for v in cofacet_vertices
        v == new_vertex && continue
        # Even though this looks like a tight loop, v changes way more often than us, so
        # this is the faster order of indexing by new_vertex and v.
        d = dist(rips, new_vertex, v)
        if ismissing(d) || d > threshold(rips)
            return nothing
        else
            _d::T = d
            diameter = ifelse(_d > diameter, _d, diameter)
        end
    end
    return S(I(sign) * index(cofacet_vertices), diameter)
end

@inline @propagate_inbounds function unsafe_cofacet(
    ::Type{S},
    rips::AbstractRipsFiltration{I},
    sx,
    cofacet_vertices,
    _,
    sign,
    new_edges,
) where {I, S<:Simplex}
    new_diam = birth(sx)
    for e in new_edges
        e > threshold(rips) && return nothing
        new_diam = ifelse(e > new_diam, e, new_diam)
    end
    new_index = index(cofacet_vertices)
    return S(I(sign) * new_index, new_diam)
end

function edges(rips::AbstractRipsFiltration)
    if issparse(dist(rips))
        _sparse_edges(rips)
    else
        _full_edges(rips)
    end
end

function _full_edges(rips::AbstractRipsFiltration)
    result = edge_type(rips)[]
    @inbounds for u in 1:size(dist(rips), 1), v in u+1:size(dist(rips), 1)
        sx = unsafe_simplex(rips, Val(1), (v, u), 1)
        !isnothing(sx) && push!(result, sx)
    end
    return result
end

function _sparse_edges(rips::AbstractRipsFiltration)
    result = edge_type(rips)[]
    rows = rowvals(rips.dist)
    vals = nonzeros(rips.dist)
    for u in 1:size(rips.dist, 1)
        for i in nzrange(dist(rips), u)
            v = rows[i]
            if v > u
                sx = unsafe_simplex(rips, Val(1), (v, u), 1)
                !isnothing(sx) && push!(result, sx)
            end
        end
    end
    return result
end

"""
    to_matrix(points)

Convert collection of d-dimensional points to d×n matrix
"""
function to_matrix(points)
    dim = length(points[1])
    T = eltype(points[1])
    n = length(points)
    result = zeros(T, (dim, n))

    for (i, p) in enumerate(points)
        length(p) == dim || throw(ArgumentError("points must have the same length"))
        result[:, i] .= p
    end

    return result
end

"""
    distances(metric, points)

Return distance matrix calculated from `points` with `metric`.
"""
function distances(metric, points)
    points_mat = to_matrix(points)
    dists = pairwise(metric, points_mat, dims=2)
    return dists
end

"""
    radius(dists)
    radius(points[, metric=Euclidean()])

Calculate the radius of the space. This is used for default `thresholds`.
"""
function radius(dists::AbstractMatrix{T}) where T
    return minimum(maximum(abs, dists[:, i]) for i in 1:size(dists, 1))
end
function radius(points, metric=Euclidean())
    radius = Inf
    for p in points
        p_max = 0.0
        for q in points
            p == q && continue
            p_max = max(p_max, metric(SVector(p), SVector(q)))
        end
        radius = min(p_max, radius)
    end
    return radius
end


"""
    Rips{I, T} <: AbstractRipsFiltration{I, T}

This type represents a filtration of Vietoris-Rips complexes.
Diagonal items are treated as vertex birth times.

Threshold defaults to radius of input space. When using low `threshold`s, consider using
[`SparseRips`](@ref) instead.

# Constructors

* `Rips(distance_matrix; threshold=nothing)`
* `Rips(points; metric=Euclidean(), threshold=nothing)`
* `Rips{I}(args...)`: `I` sets the size of integer used to represent simplices.
"""
struct Rips{I, T, A<:AbstractMatrix{T}} <: AbstractRipsFiltration{I, T}
    dist::A
    threshold::T
end

function Rips{I}(
    dist::AbstractMatrix{T}; threshold=nothing
) where {I, T}
    issymmetric(dist) ||
        throw(ArgumentError("`dist` must be symmetric"))
    !issparse(dist) ||
        throw(ArgumentError("`dist` is sparse. Use `SparseRips` instead"))
    thresh = isnothing(threshold) ? radius(dist) : T(threshold)
    return Rips{I, T, typeof(dist)}(dist, thresh)
end
function Rips{I}(points::AbstractVector; metric=Euclidean(), kwargs...) where I
    return Rips{I}(distances(metric, points); kwargs...)
end
function Rips(dist; kwargs...)
    return Rips{Int}(dist; kwargs...)
end

@propagate_inbounds function dist(rips::Rips{<:Any, T}, i, j) where T
    return rips.dist[i, j]
end
dist(rips::Rips) = rips.dist

threshold(rips::Rips) = rips.threshold
birth(rips::Rips, i) = rips.dist[i, i]

# sparse rips filtration ================================================================= #
"""
    SparseRips{I, T} <: AbstractRipsFiltration{T, Simplex}

This type represents a sparse filtration of Vietoris-Rips complexes.
The distance matrix will be converted to a sparse matrix with all values greater than
threshold deleted. Off-diagonal zeros in the matrix are treated as `missing`. Diagonal items
are treated as vertex birth times.

# Constructor

* `SparseRips{I}(distance_matrix; threshold=nothing)`
* `SparseRips(distance_matrix; threshold=nothing)`: `I` sets the integer size used to
  represent simplices.
"""
struct SparseRips{I, T} <: AbstractRipsFiltration{I, T}
    dist::SparseMatrixCSC{T, Int}
    threshold::T
end

function SparseRips{I}(
    dist::AbstractMatrix{T}; threshold=nothing
) where {I, T}
    issymmetric(dist) || throw(ArgumentError("`dist` must be symmetric"))
    if isnothing(threshold)
        threshold = issparse(dist) ? maximum(dist) : radius(dist)
    end
    dists = SparseArrays.fkeep!(SparseMatrixCSC(dist), (_, _, v) -> v ≤ threshold)

    return SparseRips{I, T}(dists, threshold)
end
function SparseRips(dist; threshold=nothing)
    return SparseRips{Int}(dist; threshold=threshold)
end
function SparseRips{I}(points::AbstractVector; metric=Euclidean(), kwargs...) where I
    return SparseRips{I}(distances(metric, points); kwargs...)
end

@propagate_inbounds function dist(rips::SparseRips{<:Any, T}, i, j) where T
    res = rips.dist[i, j]
    return ifelse(i == j, res, ifelse(iszero(res), missing, res))
end
dist(rips::SparseRips) = rips.dist

threshold(rips::SparseRips) = rips.threshold
birth(rips::SparseRips, i) = rips.dist[i, i]
simplex_type(::SparseRips{I, T}, dim) where {I, T} = Simplex{dim, T, I}

# This is the coboundary used when distance matrix in AbstractRipsFiltration is sparse.
function coboundary(
    rips::AbstractRipsFiltration, sx::AbstractSimplex, ::Val{A}=Val(true)
) where A
    if dist(rips) isa SparseMatrixCSC
        return SparseCoboundary{A}(rips, sx)
    else
        return Coboundary{A}(rips, sx)
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
        colptr = dist(filtration).colptr
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
    rowval = dist(it.filtration).rowval
    nzval = dist(it.filtration).nzval
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
