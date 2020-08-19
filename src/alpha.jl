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
        A[1:end-1, 1:end-1] .= 2 * P' * P

        b = vec(sum(abs2, P, dims=1))
        push!(b, 1.0)

        fact = lu(A, check=false)
        if issuccess(fact)
            x = fact \ b
            bary_coords = x[1:end-1]

            center = P * bary_coords
            radius = P[:, 1] - center
            return center, sum(abs2, radius)
        else
            return fill(NaN, length(pts[1])), Inf
        end
    end
end

# TODO: could be faster
function _build_dims!(dicts, triangulation, points, ::Val{D}, progress) where D
    if progress
        progbar = Progress(size(triangulation, 2), desc="Collecting $D-simplcies... ")
    end
    for face in eachcol(triangulation)
        for σ in IterTools.subsets(face, Val(D+1))
            @assert issorted(σ, rev=true)
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
            for i in 1:D+1
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
        progress && next!(progbar)
    end
end

function _fix_dim!(dicts, threshold, ::Val{D}, progress) where D
    for (idx, birth) in dicts[D + 1]
        σ = Tuple(_vertices(idx, Val(D + 1)))
        σ_idx = index(σ)
        if D > 1
            for i in 1:D+1
                τ = TupleTools.deleteat(σ, i)
                τ_idx = index(τ)
                dicts[D][τ_idx] = min(dicts[D][τ_idx], birth)
            end
        end
        corrected_birth = 2 * √dicts[D + 1][σ_idx]
        if corrected_birth > threshold
            #delete!(dicts[D + 1], σ_idx)
            dicts[D + 1][σ_idx] = corrected_birth
        else
            dicts[D + 1][σ_idx] = corrected_birth
        end
    end
end

"""
    alpha_simplices(points[, progress])

Collect all simplices and their birth times in alpha filtration.

Based on https://github.com/scikit-tda/cechmate/blob/master/cechmate/filtrations/alpha.py
"""
function alpha_simplices(points, threshold, progress, ::Type{I}) where I
    progress && printstyled("Building triangulation... ", color=:green)
    triangulation = sort!(I.(delaunay(to_matrix(points))), dims=1, rev=true)
    progress && printstyled("done.\n", color=:green)

    largest_face = tuple(maximum(eachcol(triangulation))...)
    index_overflow_check(largest_face)

    dim = length(points[1])
    dicts = [Dict{I, Float64}() for _ in 0:dim]

    for d in dim:-1:1
        _build_dims!(dicts, triangulation, points, Val(d), progress)
    end
    for i in 1:length(points)
        dicts[1][i] = 0.0
    end
    if progress
        progbar = Progress(dim, desc="Fixing birth times...     ")
    end
    for d in dim:-1:1
        _fix_dim!(dicts, threshold, Val(d), progress)
        progress && next!(progbar)
    end

    return dicts
end

"""
    Alpha{I, P<:SVector} <: AbstractFiltration{I, Float64}

`Alpha` filtrations are filtrations of the Delaunay complex.

They have much fewer simplices than `Rips`, so it's a good for large, low-dimensional
datasets.  What "low-dimensional" means depends on the data, but this is definitely a good
choice for 3D or lower. For high dimensional data, Delaunay complex construction may take a
long time.

# Note

Unlike most implementations, this one uses circumdiameters instead of circumradii. This
makes the scale of the results comparable to `Rips`. If you need radius based values, divide
your data or the resulting interval endpoints by 2.

"""
struct Alpha{I, P<:SVector} <: AbstractCustomFiltration{I, Float64}
    dicts::Vector{Dict{I, Float64}}
    adj::SparseMatrixCSC{Bool, Int}
    threshold::Float64
    points::Vector{P}
end
function Alpha{I}(points; threshold=nothing, progress=false) where I
    pts = SVector.(points)
    threshold = isnothing(threshold) ? radius(pts) : threshold
    dicts = alpha_simplices(pts, threshold, progress, I)
    adj = adjacency_matrix(dicts)
    return Alpha{I, eltype(pts)}(dicts, adj, threshold, pts)
end
function Alpha(points; kwargs...)
    return Alpha{Int}(points; kwargs...)
end

dist(alpha::Alpha) = alpha.adj
simplex_dicts(alpha::Alpha) = alpha.dicts
threshold(alpha::Alpha) = alpha.threshold

function alpha_skelly(alpha::Alpha)
    result = similar(dist(alpha), Float64)
    for e in edges(alpha)
        u, v = vertices(e)
        result[u, v] = result[v, u] = birth(e)
    end
    result
end

function alpha_distmat(points::Vector{SVector{D, Float64}}) where D
    is = Int[]
    js = Int[]
    vs = Float64[]
    faces = reinterpret(SVector{D + 1, Int32}, delaunay(to_matrix(points)))
    ncorrected=nnotcorrected=0
    corr = falses(length(faces))
    for (kk, face) in enumerate(faces)
        for i in 1:D+1, j in i+1:D+1
            others = deleteat(deleteat(face, j), i)
            p = points[face[i]]
            q = points[face[j]]
            @assert face[i] ∉ others
            @assert face[j] ∉ others

            c, r2 = circumcenter_radius2(SVector(p, q))

            if any(sum.(abs2, Ref(c) .- points[others]) .< r2)
                # When edge is close to another vertex in a facet, we can't use its radius
                # as the birth time.
                selection = sum.(abs2, Ref(c) .- points[others]) .< r2
                vec = points[others[selection]]
                append!(vec, (p, q))
                corr[kk] = true

                _, weight = circumcenter_radius2(vec)
            else
                weight = r2
            end

            append!(is, (face[i], face[j]))
            append!(js, (face[j], face[i]))
            append!(vs, (weight, weight))
        end
    end
    return sparse(is, js, sqrt.(vs), length(points), length(points), min)
end
