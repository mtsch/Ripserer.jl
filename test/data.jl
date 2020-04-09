"""
    rand_dist_matrix(n, [sparse])

Construct a random distance matrix with `n` rows and columns.
"""
function rand_dist_matrix(n)
    A = rand(n, n)
    A .+= A'
    A -= Diagonal(A)
    A
end
function rand_dist_matrix(n, sparse)
    A = sprand(n, n, sparse/2)
    A .+= A'
    A -= Diagonal(A)
    dropzeros!(A)
    A
end

"""
    torus(n)

Construct a torus distance matrix with `n` points. The points are equidistant.
"""
function torus(n)
    n = floor(Int, sqrt(n))
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

# Distances on an icosahedron graph with edge length 1.
icosahedron = Float64[0 1 2 2 1 2 1 1 2 2 1 3;
                      1 0 3 2 1 1 2 1 2 1 2 2;
                      2 3 0 1 2 2 1 2 1 2 1 1;
                      2 2 1 0 3 2 1 1 2 1 2 1;
                      1 1 2 3 0 1 2 2 1 2 1 2;
                      2 1 2 2 1 0 3 2 1 1 2 1;
                      1 2 1 1 2 3 0 1 2 2 1 2;
                      1 1 2 1 2 2 1 0 3 1 2 2;
                      2 2 1 2 1 1 2 3 0 2 1 1;
                      2 1 2 1 2 1 2 1 2 0 3 1;
                      1 2 1 2 1 2 1 2 1 3 0 2;
                      3 2 1 1 2 1 2 2 1 1 2 0]

# Distances on a cycle graph with edge length 1.
cycle = [0 1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2 1;
         1 0 1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2;
         2 1 0 1 2 3 4 5 6 7 8 9 8 7 6 5 4 3;
         3 2 1 0 1 2 3 4 5 6 7 8 9 8 7 6 5 4;
         4 3 2 1 0 1 2 3 4 5 6 7 8 9 8 7 6 5;
         5 4 3 2 1 0 1 2 3 4 5 6 7 8 9 8 7 6;
         6 5 4 3 2 1 0 1 2 3 4 5 6 7 8 9 8 7;
         7 6 5 4 3 2 1 0 1 2 3 4 5 6 7 8 9 8;
         8 7 6 5 4 3 2 1 0 1 2 3 4 5 6 7 8 9;
         9 8 7 6 5 4 3 2 1 0 1 2 3 4 5 6 7 8;
         8 9 8 7 6 5 4 3 2 1 0 1 2 3 4 5 6 7;
         7 8 9 8 7 6 5 4 3 2 1 0 1 2 3 4 5 6;
         6 7 8 9 8 7 6 5 4 3 2 1 0 1 2 3 4 5;
         5 6 7 8 9 8 7 6 5 4 3 2 1 0 1 2 3 4;
         4 5 6 7 8 9 8 7 6 5 4 3 2 1 0 1 2 3;
         3 4 5 6 7 8 9 8 7 6 5 4 3 2 1 0 1 2;
         2 3 4 5 6 7 8 9 8 7 6 5 4 3 2 1 0 1;
         1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2 1 0]

# taken from ripser/examples
projective_plane = [0 1 1 1 1 1 1 1 1 2 2 2 2;
                    1 0 2 2 2 1 2 1 2 1 2 2 2;
                    1 2 0 2 2 2 1 2 1 1 2 2 2;
                    1 2 2 0 2 1 2 2 1 2 2 2 1;
                    1 2 2 2 0 2 1 1 2 2 2 2 1;
                    1 1 2 1 2 0 2 2 2 1 1 2 1;
                    1 2 1 2 1 2 0 2 2 1 1 2 1;
                    1 1 2 2 1 2 2 0 2 1 2 1 1;
                    1 2 1 1 2 2 2 2 0 1 2 1 1;
                    2 1 1 2 2 1 1 1 1 0 1 1 2;
                    2 2 2 2 2 1 1 2 2 1 0 2 1;
                    2 2 2 2 2 2 2 1 1 1 2 0 1;
                    2 2 2 1 1 1 1 1 1 2 1 1 0]
