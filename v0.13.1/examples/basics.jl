# # Basics

# In this example we will present the usage of Ripserer. We start by loading some packages.

using Ripserer
using Plots
using Random; Random.seed!(1337); gr(); nothing # hide

# ## Basic Usage

# Let's start with a basic example, points randomly sampled from a noisy circle.
# We start by defining our sampling function.

function noisy_circle(n; r=1, noise=0.1)
    points = NTuple{2, Float64}[]
    for _ in 1:n
        θ = 2π * rand()
        push!(points, (r*sin(θ) + noise*rand(), r*cos(θ) + noise*rand()))
    end
    points
end

# Next, we sample 100 points from the circle.

circ_100 = noisy_circle(100)
scatter(circ_100, aspect_ratio=1, legend=false, title="Noisy Circle")

# To compute the persistent homology, simply run the following. The `dim_max` argument sets
# the maximum dimension persistent homology is computed in.

result_circ = ripserer(circ_100, dim_max=3)

# !!! warning "Warning"
#     Computing Vietoris-Rips persistent homology in high dimensions for large numbers of
#     points is computationally expensive and requires a large amount of memory. Be careful
#     or you **will** run out of memory. On an ordinary computer, you can expect to compute
#     one-dimensional persistent homology for datasets of a few thousand points and higher
#     (2-3) dimensional persistent homology for datasets of a few hundred points. This, of
#     course, depends on the data set itself.

# The result can be plotted as a persistence diagram.

plot(result_circ)

# Or as a barcode.

barcode(result_circ)

# ``H_1``, ``H_2`` and ``H_3`` in this plot are hard to see, because we have too many
# ``H_0`` bars. We can plot only some of the diagrams.

# !!! note "Note"
#     `result` is just an array of persistence diagrams, so the zero-dimensional diagram is
#     found at index 1.

barcode(result_circ[2:end], linewidth=2)

# We can plot a single diagram in the same manner.

barcode(result_circ[3], linewidth=3)

# ## Distance Matrix Inputs

# In the previous example, we got our result by passing a collection of points to
# `ripserer`.  Under the hood, the algorithm actually works with distance matrices. Let's
# define a distance matrix of the shortest paths on a [regular
# icosahedron](https://en.wikipedia.org/wiki/Regular_icosahedron) graph.

# ```@raw html
# <img src="https://upload.wikimedia.org/wikipedia/commons/8/83/Icosahedron_graph.svg" height="250" width="250">
# ```

icosahedron = [0 1 2 2 1 2 1 1 2 2 1 3;
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
nothing; # hide

# To compute the persistent homology, simply feed the distance matrix to `ripserer`.

result_icosa = ripserer(icosahedron, dim_max=2)

# Because an icosahedron is topologically equivalent to a sphere, we got a single class in
# the second dimension.

result_icosa[3]

# ## Thresholding

# In our next example, we will show how to use thresholding to speed up computation. We
# start by defining a sampling function that generates ``n`` points from the square
# ``[-4,4]\times[-4,4]`` with a circular hole of radius 1 in the middle.

function cutout(n)
    points = NTuple{2, Float64}[]
    while length(points) < n
        x, y = (8rand() - 4, 8rand() - 4)
        if x^2 + y^2 > 1
            push!(points, (x, y))
        end
    end
    points
end

# We sample 2000 points from this space.

cutout_2000 = cutout(2000)
scatter(cutout_2000, markersize=1, aspect_ratio=1, legend=false, title="Cutout")

# We calculate the persistent homology and time the calculation.

@time result_cut = ripserer(cutout_2000)
nothing # hide

#

plot(result_cut)

# Notice that while there are many 1-dimensional classes, one of them stands out. This class
# represents the hole in the middle of our square. We can extract this interval by doing the
# following.

most_persistent = sort(result_cut[2], by=persistence)[end]

# Notice the death time of this interval is around 1.83 and that no intervals occur after
# that time. This means that we could stop computing when we reach this time and the result
# should not change. Let's try it out!

@time result_cut_thresh_2 = ripserer(cutout_2000, threshold=2)
nothing; # hide

#

plot(result_cut_thresh_2, title="Persistence Diagram, threshold=2")

# Indeed, the result is exactly the same, but it took less than a third of the time to
# compute.

result_cut_thresh_2 == result_cut

# If we pick a threshold that is too low, we still detect the interval, but its death time
# becomes infinite.

@time result_cut_thresh_1 = ripserer(cutout_2000, threshold=1)
nothing; # hide

#

plot(result_cut_thresh_1, title="Persistence Diagram, threshold=1")
