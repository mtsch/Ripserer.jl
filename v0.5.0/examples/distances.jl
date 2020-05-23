# # Distances

# In this example we will demonstrate computing distances between persistence diagrams.

using Ripserer
using Plots
using Random; Random.seed!(1337); gr(); nothing # hide

# We will again look at the persistent homology of a square with a round hole in the center.

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
nothing # hide

# In this example, we're interested in the first persistent homology group of two samples of
# this space, one with 1000 points and another with 500 points.

cutout_1 = cutout(1000)
scatter(cutout_1,
        color=2,
        aspect_ratio=1,
        label="",
        title="cutout_1")

#

cutout_2 = cutout(500)
scatter(cutout_2,
        color=3,
        aspect_ratio=1,
        label="",
        title="cutout_2")

# We compute the first persistence diagrams.

_, result_1 = ripserer(cutout_1)
_, result_2 = ripserer(cutout_2)

plot(plot(result_1, title="cutout_1", markercolor=2, xlim=(0, 2.2), ylim=(0, 2.2)),
     plot(result_2, title="cutout_2", markercolor=3, xlim=(0, 2.2), ylim=(0, 2.2)))

# We notice the diagrams are similar. We can formally describe this similarity with the
# Bottleneck and Wasserstein distances. First, let's look at their definitions.

# The ``q``-th Wasserstein distance is defined as

# ```math
# W_q(X,Y)=\left[\inf_{\eta:X\rightarrow Y}\sum_{x\in X}||x-\eta(x)||_\infty^q\right],
# ```

# where ``X`` and ``Y`` are the persistence diagrams and ``\eta`` is a perfect matching
# between the intervals. We also include the points on the diagonals of ``X`` and ``Y`` to
# ensure a perfect matching always exists.

# The bottleneck distance ``W_\infty`` is defined in a similar manner, as

# ```math
# W_\infty(X, Y) = \inf_{\eta:X\rightarrow Y} \sup_{x\in X} ||x-\eta(x)||_\infty.
# ```

# In other words, the value of the Wasserstein distance between to diagrams is equal to the
# ``q``-norm of the of the solution of the [assignment
# problem](https://en.wikipedia.org/wiki/Assignment_problem). The value of the Bottleneck
# distance is equal to the maximum weight in the solution to the minimum weight perfect
# matching problem.

# Now, let's look at what these distances say about the diagrams we have computed earlier.

# We can construct a matching by using the `matching` function.

match_bottle = matching(Bottleneck(), result_1, result_2)

#

match_wasser = matching(Wasserstein(), result_1, result_2)

# We can access the distance as follows.

distance(match_bottle)

# Or if we don't care about the matching, we can compute the distance directly.

distance(Wasserstein(), result_1, result_2)

# We plot the matchings using the `plot` function.

plot(match_bottle)

#

plot(match_wasser)

# In the case of the bottleneck distance, we can show the full matching by supplying the
# `bottleneck=false` parameter to `plot`.

plot(match_bottle, bottleneck=false)
