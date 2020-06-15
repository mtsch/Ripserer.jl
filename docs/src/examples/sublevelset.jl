# # Sublevel Set Persistent Homology

# In this example, we will demonstrate computing sublevel set persistent homology of time
# series and image data. We will need the following packages.

using Ripserer
using Images
using Plots
gr() # hide
using Random
Random.seed!(1337)

# We will use two data sets in this example. The first will be a curve:

n = 1000
x = range(0, 1, length=n)
curve = sin.(2π * 5x) .* x

curve_plot = plot(curve, legend=false, title="Curve")

# The other will be [the Event Horizon Telescope picture of a black
# hole](https://en.wikipedia.org/wiki/File:Black_hole_-_Messier_87_crop_max_res.jpg). We
# will use a small, 240×240 pixel version of the image. Ripserer should have no problems
# with processing larger images, but this will work well enough for the purpose of this
# tutorial.

blackhole_image = load(
    joinpath(@__DIR__, "../assets/data/240px-Black_hole_-_Messier_87_crop_max_res.jpg")
)
blackhole_plot = plot(blackhole_image, title="Black Hole")

# In order to use the image with Ripserer, we have to convert it to grayscale.

blackhole = Gray.(blackhole_image)
nothing # hide

# ## One-dimensional Case

# Sublevel set persistent homology provides a stable descriptions of the critical points of
# a function. The zeroth persistent homology group ``H_0`` corresponds to the local minima
# of a function. To compute this with Ripserer, we use cubical persistent homology as
# follows. Note that there is not information in ``H_1``, since the function is
# one-dimensional.

result, _ = ripserer(Cubical(curve))
plot(curve_plot, plot(result))

# To see which minimum each interval corresponds to, we compute representatives.

result, _ = ripserer(Cubical(curve), reps=true)

# The infinite interval's representative will include the whole function. To plot a
# representative, simply pass the interval along with the data to `plot`. An alternative
# way to plot the same thing is to get the representative of the interval first and then
# plot that in the same way. Let's use both ways in the following example.

infinite_interval = only(filter(!isfinite, result))
plt = plot(representative(infinite_interval), curve,
           legend=false,
           title="Representatives",
           seriestype=:path)

for interval in filter(isfinite, result)
    plot!(plt, interval, curve, seriestype=:path)
end

# To get the locations of the minima, extract the critical simplices from intervals. As
# simplices act like collections of vertex indices, we can use `ònly` to get them.

x_min = [first(birth_simplex(int)) for int in result]

scatter!(plt, x_min, curve[x_min], color=1:6, markershape=:star)
plot(plt, plot(result, markercolor=1:6, markeralpha=1))


# Note that each interval's birth is equal to the value of the corresponding local minimum
# and its death is equal to the higher of the two adjacient maxima. An intuitive way of
# thinking about the result is imagining you pour water in the curve. Water is collected in
# a valley and once it reaches a local maximum, it starts pouring in the adjacient valley.

# ## Two-dimensional Case

# Now let's do a similar thing for a 2d example. There is nothing stopping us from going
# into higher dimensions, but we will skip those. The principles are exactly the same.
# Instead of looking for local minima, let's look for local maxima. To do that, we have to
# invert the image.

result = ripserer(Cubical(1 .- blackhole))

# !!! tip "Overflow"
#     Ripserer stores vertices in large integers. In some cases, `Int64` is
#     not large enough to represent them. This is particularly likely to happen with
#     `Cubical`. When it does, `ripserer` will throw an `OverflowError`. If that happens,
#     try constructing the `Cubical` filtration with `Cubical{Int128}(data)`.

plot(blackhole_plot, plot(result))

# We notice there are quite a lot of intervals along the diagonal. These correspond to local
# geometry of the image, so we are not interested in them right now. To filter them out, we
# set a cutoff. All intervals with persistence lower than cutoff will be ignored.

result = ripserer(Cubical(1 .- blackhole), cutoff=0.1)
plot(blackhole_plot, plot(result))

# Like earlier, we can show the local extrema in the image. We will show a different way to
# plot them. We use the `threshold` argument with `plot`, which only keeps parts of the
# representative with a diameter equal to or lower than `threshold`. If we needed a strict
# `<`, we could use `threshold_strict`.

result = ripserer(Cubical(1 .- blackhole), cutoff=0.1, reps=true)
plt = plot(blackhole_image, title="Black Hole")
for interval in result[1]
    plot!(plt, interval, blackhole, threshold=birth(interval), label=string(interval))
end
plot(plt, plot(result[1], markercolor=2:3, markeralpha=1))

# Some maxima have multiple values because more than one pixel in the image has the same
# value.

# Unlike with the previous example, we now also have access to ``H_1``, which corresponds to
# the cycles in the image. Let's try to plot the representative of ``H_1``.

plt = plot(blackhole_image, title="Black Hole")
plot!(plt, only(result[2]), blackhole, label="H₁ cocycle", color=1)

# Notice that the result is not actually a cycle, but rather a collection of pixels that
# would destroy the cycle if removed. The reason is that Ripserer computes persistent
# _co_homology by default. The persistence diagrams of persistent homology and persistent
# cohomology are the same, but persistent cohomology is much more efficient to compute. The
# representatives it finds, however tend to not be as informative. Keep this in mind when
# trying persistent homology out for larger datasets; it can take a very very long
# time. This case is quite small and computing persistent homology should pose no problem.
# We compute persistent homology with the argument `cohomology=false`.

# !!! warning "Homology"
#     Ripserer currently can't compute infinite intervals in dimensions
#     higher than zero with persistent homology.

result = ripserer(Cubical(1 .- blackhole), cutoff=0.1, reps=true, cohomology=false)
plot!(plt, only(result[2]), blackhole, label="H₁ cycle", color=3)

# We have successfuly found the hole in a black hole.
