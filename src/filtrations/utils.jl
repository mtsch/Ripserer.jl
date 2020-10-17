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
    distances(points, metric=Euclidean(1e-12))

Return distance matrix calculated from `points` with `metric`.
"""
function distances(points, metric=Euclidean(1e-12))
    points_mat = to_matrix(points)
    dists = pairwise(metric, points_mat; dims=2)
    return dists
end

"""
    radius(dists)
    radius(points[, metric=Euclidean(1e-12)])

Calculate the radius of the space. This is used for default `thresholds`.
"""
function radius(dists::AbstractMatrix)
    return minimum(maximum(abs, dists[:, i]) for i in 1:size(dists, 1))
end
function radius(dists::SparseMatrixCSC)
    return maximum(dists)
end
function radius(points, metric=Euclidean(1e-12))
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
