# # Sublevel Image Filtrations


# In this example, we will demonstrate sublevel persistent homology on images. We will use the
# [picture of a black hole](https://en.wikipedia.org/wiki/File:Black_hole_-_Messier_87_crop_max_res.jpg).
# We will use the 1024×1024 pixel image.

# Start by loading the required packages and the image.

using Ripserer
using Plots
using FileIO
gr(); nothing # hide
blackhole = load(joinpath(@__DIR__, "../assets/blackhole1024px.jpg"))

plot(blackhole)

# Extract the red channel from the image.
blackhole_red = getfield.(blackhole, :r)

# And plot it as a heatmap.
heatmap(blackhole_red, aspect_ratio=1)

# Cumpute the persistent homology.

# !!! note "Note"
#     Unlike with Rips filtrations, we can use lare data sets here because the number of
#     simplices in a cubical filtration is linear to the size of the input. We still want to
#     be careful with calculating representatives because calculating representatives for
#     thousands of persistence intervals can take a while.

result_min = ripserer(CubicalFiltration(blackhole_red))

#

plot(result_min, markersize=1)

# We notice there is a lot of noise along the diagonal, but three intervals stand out. We
# remove the noise by supplying the `cutoff` argument to `ripserer`. This removes all
# intervals with persistence lower than cutoff.

result_min = ripserer(CubicalFiltration(blackhole_red), cutoff=0.1)

#

plot(result_min)

# Now that we have a smaller number of intervals, we can compute the representatives.
result_min = ripserer(CubicalFiltration(blackhole_red), cutoff=0.1, representatives=true)

# We separate the finite and infinite intervals. The finite one corresponds to the local
# minimum in the middle of the image and the infinite one corresponds to the global minimum.
finite_interval = only(filter(isfinite, result_min[1]))
infinite_interval = only(filter(!isfinite, result_min[1]))

# Plot the location of the minima on the image. Note that we have to use the `threshold`
# argument to only plot the simplices with diameter lower than or equal to the birth of the
# interval. For the global one, we will only color the pixels that have a value that equals
# to the interval's birth time.
heatmap(blackhole_red, aspect_ratio=1)
plot!(finite_interval, blackhole_red,
      markershape=:star,
      threshold=birth(finite_interval),
      label="local minimum")

plot!(infinite_interval, blackhole_red,
      markersize=0,
      threshold=birth(infinite_interval),
      label="global minima")

# To plot the region of the local minimum, we skip the threshold and choose a small marker
# size. We skip the global minimum region because it encompasses the whole image.
heatmap(blackhole_red, aspect_ratio=1)
plot!(finite_interval, blackhole_red,
      markersize=0,
      alpha=0.1,
      color=1,
      label="local minimum region")

# We can also plot the representative of ``H_1``, but keep in mind we have computed
# persistent cohomology. The representative gives us a chain that kills the first homology
# group.

heatmap(blackhole_red, aspect_ratio=1)
scatter!(result_min[2][1], blackhole_red,
         color=1,
         label="H₁ killer",
         markersize=0,
         linewidth=4)

# Let's do the same thing for maxima. We do this by simply negating the image.

result_max = ripserer(CubicalFiltration(-blackhole_red), cutoff=0.1, representatives=true)

heatmap(blackhole_red, aspect_ratio=1)
scatter!(result_max[2][1], blackhole_red,
         color=1,
         label="H₁ killer",
         markersize=0,
         linewidth=4)
plot!(result_max[1][1], blackhole_red,
      markersize=0,
      color=3,
      alpha=0.1,
      threshold=birth(result_max[1][1]),
      label="global maxima")

# Note that the noise in global maxima is due to jpg compression.
