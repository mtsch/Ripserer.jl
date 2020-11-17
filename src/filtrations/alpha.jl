"""
    circumcenter_radius2(pts)

Calculate circumcenter and circumradius squared from points.

Based on https://github.com/hirani/pydec/blob/master/pydec/math/circumcenter.py and
http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
"""
# TODO: try this with SMatrix?
function circumcenter_radius2(pts)
    if length(pts) == 1
        return pts[1], 0.0
    elseif length(pts) == 2
        return 0.5 * sum(pts), 0.25 * sum(abs2, pts[1] - pts[2])
    else
        P = to_matrix(pts)
        n = length(pts)

        A = ones(n + 1, n + 1)
        A[end, end] = 0.0
        A[1:(end - 1), 1:(end - 1)] .= 2 * P' * P

        b = vec(sum(abs2, P; dims=1))
        push!(b, 1.0)

        fact = lu(A; check=false)
        if issuccess(fact)
            x = fact \ b
            bary_coords = x[1:(end - 1)]

            center = P * bary_coords
            radius = P[:, 1] - center
            return center, sum(abs2, radius)
        else
            return fill(NaN, length(pts[1])), Inf
        end
    end
end

# TODO: could be faster, takes a long time to compile
function _build_dims!(dicts, triangulation, points, ::Val{D}, verbose) where {D}
    if verbose
        progbar = Progress(size(triangulation, 2); desc="Collecting $D-simplcies... ")
    end
    for face in eachcol(triangulation)
        for σ in IterTools.subsets(face, Val(D + 1))
            @assert issorted(σ; rev=true)
            σ_idx = index(σ)
            if !haskey(dicts[D + 1], σ_idx)
                _, σ_r2 = circumcenter_radius2(points[SVector(σ)])
                if isfinite(σ_r2)
                    dicts[D + 1][σ_idx] = σ_r2
                else
                    continue
                end
            end
            # Propagate birth time to facets.
            σ_r2 = dicts[D + 1][σ_idx]
            for i in 1:(D + 1)
                τ = TupleTools.deleteat(σ, i)
                τ_idx = index(τ)
                if haskey(dicts[D], τ_idx)
                    dicts[D][τ_idx] = min(σ_r2, dicts[D][τ_idx])
                elseif length(τ) > 1
                    τ_c, τ_r2 = circumcenter_radius2(points[SVector(τ)])
                    if sum(abs2, τ_c - points[σ[i]]) < τ_r2
                        dicts[D][τ_idx] = σ_r2
                    end
                end
            end
        end
        verbose && next!(progbar)
    end
end

function _fix_dim!(dicts, threshold, ::Val{D}, verbose) where {D}
    for (idx, birth) in dicts[D + 1]
        σ = Tuple(_vertices(idx, Val(D + 1)))
        σ_idx = index(σ)
        if D > 1
            for i in 1:(D + 1)
                τ = TupleTools.deleteat(σ, i)
                τ_idx = index(τ)
                dicts[D][τ_idx] = min(dicts[D][τ_idx], birth)
            end
        end
        corrected_birth = 2 * √dicts[D + 1][σ_idx]
        dicts[D + 1][σ_idx] = corrected_birth
    end
end

"""
    alpha_simplices(points[, verbose])

Collect all simplices and their birth times in alpha filtration.

Based on https://github.com/scikit-tda/cechmate/blob/master/cechmate/filtrations/alpha.py
"""
function alpha_simplices(points, threshold, verbose, ::Type{I}) where {I}
    @prog_print verbose "Building triangulation... "
    triangulation = I.(delaunay(to_matrix(points)))
    sort!.(eachcol(triangulation), rev=true, alg=InsertionSort)
    @prog_println verbose "done."

    largest_face = tuple(maximum(eachcol(triangulation))...)
    index_overflow_check(largest_face)

    dim = length(points[1])
    dicts = [Dict{I,Float64}() for _ in 0:dim]

    # Build the filtration
    for d in dim:-1:1
        _build_dims!(dicts, triangulation, points, Val(d), verbose)
    end
    for i in 1:length(points)
        dicts[1][i] = 0.0
    end
    if verbose
        progbar = Progress(dim; desc="Fixing birth times...     ")
    end
    # Make sure all simplices are born after their facets and sqrt the birth times.
    for d in dim:-1:1
        _fix_dim!(dicts, threshold, Val(d), verbose)
        verbose && next!(progbar)
    end

    return dicts
end

"""
    Alpha{I, P<:SVector} <: AbstractFiltration{I, Float64}

`Alpha` filtrations are filtrations of the Delaunay complex.

They have much fewer simplices than `Rips`, so they are efficient even with large datasets,
as long as their dimensionality is low.  What "low" means depends on the data, but this is
definitely a good choice for 3D or lower. For high dimensional data, filtration construction
may take a long time.

!!! note
    Unlike most implementations, this one uses circumdiameters instead of
    circumradii. This makes the scale of the results comparable to `Rips`. If you need
    radius based values, divide your data or the resulting interval endpoints by 2.

!!! warning
    This filtration uses [MiniQhull.jl](https://github.com/gridap/MiniQhull.jl). Please see
    the installation instructions if constructions cause errors. MiniQhull currently has
    problems running on Windows. See [this
    issue](https://github.com/gridap/MiniQhull.jl/issues/5) for more info.

# Constructors

* `Alpha(points; threshold, verbose)`: `points` should be a vector of `Tuple`s, `SVector`s
  or similar.
* `Alpha{I}(args...)`: `I` sets the size of integer used to represent simplices. Try using
  `I=Int128` if construction complains about overflow.

# Reference

Edelsbrunner, H. (1993, July). The union of balls and its dual shape. [In Proceedings of the
ninth annual symposium on Computational geometry
(pp. 218-231)](https://dl.acm.org/doi/abs/10.1145/160985.161139).

# Example

```jldoctest
julia> data = [(sin(t), cos(t), (t - π)^2) for t in range(0, 2π, length=101)[1:end-1]];

julia> alpha = Alpha(data)
Alpha{Int64, Float64}(nv=100)

julia> rips = Rips(data)
Rips{Int64, Float64}(nv=100, sparse=false)

julia> length(Ripserer.edges(alpha))
197

julia> length(Ripserer.edges(rips))
3613

julia> sort(ripserer(alpha)[2], by=persistence)[end]
[0.375, 2.01) with:
 birth_simplex: Ripserer.Simplex{1,Float64,Int64}
 death_simplex: Ripserer.Simplex{2,Float64,Int64}

julia> sort(ripserer(rips)[2], by=persistence)[end]
[0.375, 2.01) with:
 birth_simplex: Ripserer.Simplex{1,Float64,Int64}
 death_simplex: Ripserer.Simplex{2,Float64,Int64}
```
"""
struct Alpha{I,P<:SVector} <: AbstractCustomFiltration{I,Float64}
    dicts::Vector{Dict{I,Float64}}
    adj::SparseMatrixCSC{Bool,Int}
    threshold::Float64
    points::Vector{P}
end
function Alpha{I}(points; threshold=nothing, verbose=false) where {I}
    pts = SVector.(points)
    threshold = isnothing(threshold) ? 2radius(pts) : threshold
    dicts = alpha_simplices(pts, threshold, verbose, I)
    adj = _adjacency_matrix(dicts)
    return Alpha{I,eltype(pts)}(dicts, adj, threshold, pts)
end
function Alpha(points; kwargs...)
    return Alpha{Int}(points; kwargs...)
end

adjacency_matrix(alpha::Alpha) = alpha.adj
simplex_dicts(alpha::Alpha) = alpha.dicts
threshold(alpha::Alpha) = alpha.threshold

struct AlphaDist{A} <: AbstractMatrix{Float64}
    alpha::A
end

Base.size(ad::AlphaDist) = (nv(ad.alpha), nv(ad.alpha))
function Base.getindex(ad::AlphaDist, i::Integer, j::Integer)
    return Euclidean()(ad.alpha.points[i], ad.alpha.points[j])
end

distance_matrix(alpha::Alpha) = AlphaDist(alpha)
