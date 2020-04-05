using Distances
using Random
using LinearAlgebra

function torus(n)
    r = range(-1, 1, length = n+1)[1:end-1]
    pts = [(i, j) for i in r for j in r]
    dist = fill(Inf, length(pts), length(pts))
    for i in 1:length(pts), j in i+1:length(pts)
        for x_offset in -2:2:2, y_offset in -2:2:2
            x1, y1 = pts[i]
            x2, y2 = pts[j] .+ (x_offset, y_offset)
            dist[i, j] = dist[j, i] = min(√((x1-x2)^2 + (y1-y2)^2), dist[i, j])
        end
    end
    for i in 1:length(pts)
        dist[i, i] = 0
    end
    dist
end

function rand_sphere(n, dim=2)
    pts = mapslices(normalize, rand(dim, n), dims=1)
    pairwise(Euclidean(), pts)
end

dist_torus = torus(10)
