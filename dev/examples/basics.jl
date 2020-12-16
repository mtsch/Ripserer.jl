# # Usage Guide

# In this example, we will present the basics of using Ripserer. We start by loading some
# packages.

using Distances
using Plots
using Ripserer
using Random # hide
Random.seed!(1337) # hide
gr() # hide
nothing # hide

# ## Using Ripserer With Point Cloud Data

# Let's start with generating some points, randomly sampled from a noisy circle.

function noisy_circle(n; r=1, noise=0.1)
    points = NTuple{2,Float64}[]
    for _ in 1:n
        θ = 2π * rand()
        push!(points, (r * sin(θ) + noise * rand(), r * cos(θ) + noise * rand()))
    end
    return points
end

circ_100 = noisy_circle(100)
scatter(circ_100; aspect_ratio=1, legend=false, title="Noisy Circle")

# !!! tip "Point-like data types"
#     Ripserer can interpret various kinds of data as point clouds. The limitation is that
#     the data set should be an `AbstractVector` with elements with the following
#     properties:
#     * all elements are collections numbers;
#     * all elements have the same length.
#     Examples of element types that work are `Tuple`s,
#     [`SVector`](https://github.com/JuliaArrays/StaticArrays.jl)s, and
#     [`Point`](https://github.com/JuliaGeometry/GeometryBasics.jl)s.

# To compute the Vietoris-Rips persistent homology of this data set, run the
# following.

ripserer(circ_100)

# You can use the `dim_max` argument to set the maximum dimension persistent homology is
# computed in.

result_rips = ripserer(circ_100; dim_max=3)

# The result can be plotted as a persistence diagram or as a barcode.

plot(result_rips)
barcode(result_rips)
plot(plot(result_rips), barcode(result_rips)) # hide

# We can also plot a single diagram or a subset of all diagrams in the same manner. Keep in
# mind that the result is just a vector of [`PersistenceDiagram`](@ref)s. The
# zero-dimensional diagram is found at index 1.

plot(result_rips[2])
barcode(result_rips[2:end]; linewidth=2)
plot(plot(result_rips[2]), barcode(result_rips[2:end]; linewidth=2)) # hide

# Plotting can be further customized using the standard attributes from
# [Plots.jl](http://docs.juliaplots.org/latest/).

plot(result_rips; markeralpha=1, markershape=:star, color=[:red, :blue, :green, :purple])

# ## Changing Filtrations

# By default, calling [`ripserer`](@ref) will compute persistent homology with the
# [`Rips`](@ref) filtration. To use a different filtration, we have two options.

# The first option is to pass the filtration constructor as the first argument. Any keyword
# arguments the filtration accepts can be passed to [`ripserer`](@ref) and it will be
# forwarded to the constructor.

ripserer(EdgeCollapsedRips, circ_100; threshold=1, dim_max=3, metric=Euclidean())

# The second option is to initialize the filtration object first and use that as an argument
# to [`ripserer`](@ref). This can be useful in cases where constructing the filtration takes
# a long time.

collapsed_rips = EdgeCollapsedRips(circ_100; threshold=1, metric=Euclidean())
ripserer(collapsed_rips; dim_max=3)

# ## Distance Matrix Inputs

# In the previous example, we got our result by passing a collection of points to
# [`ripserer`](@ref). Under the hood, [`Rips`](@ref) and [`EdgeCollapsedRips`](@ref)
# actually work with distance matrices. Let's define a distance matrix of the shortest
# paths on a [regular icosahedron](https://en.wikipedia.org/wiki/Regular_icosahedron) graph.

# ```@raw html
# <img src="https://upload.wikimedia.org/wikipedia/commons/8/83/Icosahedron_graph.svg" height="200" width="200">
# ```

icosahedron = [
    0 1 2 2 1 2 1 1 2 2 1 3
    1 0 3 2 1 1 2 1 2 1 2 2
    2 3 0 1 2 2 1 2 1 2 1 1
    2 2 1 0 3 2 1 1 2 1 2 1
    1 1 2 3 0 1 2 2 1 2 1 2
    2 1 2 2 1 0 3 2 1 1 2 1
    1 2 1 1 2 3 0 1 2 2 1 2
    1 1 2 1 2 2 1 0 3 1 2 2
    2 2 1 2 1 1 2 3 0 2 1 1
    2 1 2 1 2 1 2 1 2 0 3 1
    1 2 1 2 1 2 1 2 1 3 0 2
    3 2 1 1 2 1 2 2 1 1 2 0
]
nothing # hide

# To compute the persistent homology, simply feed the distance matrix to [`ripserer`](@ref).

result_icosa = ripserer(icosahedron; dim_max=2)

# ## Thresholding

# In our next example, we will show how to use thresholding to speed up computation. We
# start by defining a sampling function that generates ``n`` points from the square
# ``[-4,4]\times[-4,4]`` with a circular hole of radius 1 in the middle.

function cutout(n)
    points = NTuple{2,Float64}[]
    while length(points) < n
        x, y = (8rand() - 4, 8rand() - 4)
        if x^2 + y^2 > 1
            push!(points, (x, y))
        end
    end
    return points
end

# We sample 2000 points from this space.

cutout_2000 = cutout(2000)
scatter(cutout_2000; markersize=1, aspect_ratio=1, legend=false, title="Cutout")

# We calculate the persistent homology and time the calculation.

@time result_cut = ripserer(cutout_2000)
nothing # hide

#

plot(result_cut)

# Notice that while there are many 1-dimensional classes, one of them stands out. This class
# represents the hole in the middle of our square. Since the intervals are sorted by
# persistence, we know the last interval in the diagram will be the most persistent.

most_persistent = result_cut[2][end]

# Notice the death time of this interval is around 1.83 and that no intervals occur after
# that time. This means that we could stop computing when we reach this time and the result
# should not change. Let's try it out.

@time result_cut_thresh_2 = ripserer(cutout_2000; threshold=2)
nothing # hide

#

plot(result_cut_thresh_2; title="Persistence Diagram, threshold=2")

# Indeed, the result is exactly the same, but it took less than a third of the time to
# compute.

@assert result_cut_thresh_2 == result_cut # hide
result_cut_thresh_2 == result_cut

# If we pick a threshold that is too low, we still detect the interval, but its death time
# becomes infinite.

@time result_cut_thresh_1 = ripserer(cutout_2000; threshold=1)
nothing # hide

#

result_cut_thresh_1[2][end]

# ## Persistence Diagrams

# The result of a computation is returned as a vector of
# [`PersistenceDiagram`](@ref)s. Let's take a closer look at one of those.

diagram = result_cut[2]

# The diagram is a structure that acts as a vector of
# [`PersistenceInterval`](@ref)s. As such, you can use standard Julia
# functions on the diagram.

# For example, to extract the last three intervals by birth time, you can do
# something like this.

sort(diagram; by=birth, rev=true)[1:3]

# To find the [`persistence`](@ref)s of all the intervals, you can use broadcasting.

persistence.(diagram)

# Unlike regular vectors, a [`PersistenceDiagram`](@ref) has additional metadata attached
# to it. To see all metadata, use
# [`propertynames`](https://docs.julialang.org/en/v1/base/base/#Base.propertynames).

propertynames(diagram)

# You can access the properties with the dot syntax.

diagram.field

# The attributes `dim` and `threshold` are given special treatment and can be extracted
# with appropriately named functions.

dim(diagram), threshold(diagram)

# Now, let's take a closer look at one of the intervals.

interval = diagram[end]

# An interval is very similar to a tuple of two `Float64`s, but also has some metadata
# associated with it.

interval[1], interval[2]

# [`birth`](@ref), [`death`](@ref), [`persistence`](@ref), and [`midlife`](@ref) can be used
# to query commonly used values.

birth(interval), death(interval), persistence(interval), midlife(interval)

# Accessing metadata works in a similar manner as with diagrams.

propertynames(interval)

#

interval.birth_simplex

#

interval.death_simplex

# ## Simplices

# In the previous section, we saw each interval has an associated [`birth_simplex`](@ref)
# and [`death_simplex`](@ref). These values are of the type [`Simplex`](@ref). Let's take a
# closer look at simplices.

simplex = interval.death_simplex

# [`Simplex`](@ref) is an internal data structure that uses some tricks to increase
# efficiency. For example, if we were to
# [`dump`](https://docs.julialang.org/en/v1/base/io-network/#Base.dump) it, we notice the
# vertices are not actually stored in the simplex itself.

dump(simplex)

# To access the vertices, we use [`vertices`](@ref).

vertices(simplex)

# Other useful attributes a simplex has are [`index`](@ref), [`dim`](@ref), and
# [`birth`](@ref).

index(simplex), dim(simplex), birth(simplex)

# A few additional notes on simplex properties.

# * A `D`-dimensional simplex is of type `Simplex{D}` and has `D + 1` vertices.
# * [`vertices`](@ref) are always sorted in descending order.
# * [`index`](@ref) and [`dim`](@ref) can be used to uniquely identify a given simplex.
# * [`birth`](@ref) determines when a simplex is added to a filtration.

# ## Conclusion

# This concludes the basic usage of Ripserer. For more detailed information, please check
# out the [API](@ref) page, as well as other examples.
