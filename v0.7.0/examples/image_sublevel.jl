# # Sublevel Image Filtrations


# In this example, we will demonstrate sublevel persistent homology on images. We will use
# the [Horizon telescope picture of a black
# hole](https://en.wikipedia.org/wiki/File:Black_hole_-_Messier_87_crop_max_res.jpg). We
# will use the 768×768 pixel image. Computing the same result with a larger image should not
# be a problem, but the plots would make this page less responsive.

# Start by loading the required packages and the image.

using Ripserer
using Plots
using Images
gr(); nothing # hide
blackhole = load(joinpath(@__DIR__, "../assets/blackhole768px.jpg"))

plot(blackhole, size=(800, 800))

# Convert the image to a matrix of floats.
blackhole_grayscale = Float64.(Gray.(blackhole))

# And plot it as a heatmap. Note that the image is flipped because matrices start at the top
# and images start at the bottom.

heatmap(blackhole_grayscale, aspect_ratio=1, size=(800, 800))

# Compute the persistent homology.

# !!! note "Note"
#     Unlike with Rips filtrations, we can use large data sets here because the number of
#     simplices in a cubical filtration is linear to the size of the input. We still want to
#     be careful with computing representatives because collecting them for thousands of
#     persistence intervals can take a while.

result_min = ripserer(Cubical(blackhole_grayscale))

# We plot the diagram with the `persistence` argument. This plots persistence vs birth
# instead of death vs birth.

plot(result_min, persistence=true)

# We notice there is a lot of noise along the diagonal, but four intervals stand out. We
# remove the noise by supplying the `cutoff` argument to `ripserer`. This removes all
# intervals with persistence strictly smaller than cutoff.

result_min = ripserer(Cubical(blackhole_grayscale), cutoff=0.1)

#

plot(result_min, persistence=true)

# Now that we know we have a smaller number of intervals, we can compute the
# representatives.
result_min = ripserer(Cubical(blackhole_grayscale), cutoff=0.1, reps=true)

# We separate the finite and infinite intervals. The finite one corresponds to the local
# minimum in the middle of the image and the infinite one corresponds to the global minimum.
finite_interval = only(filter(isfinite, result_min[1]))
infinite_interval = only(filter(!isfinite, result_min[1]))
nothing; # hide

# Plot the location of the minima on the image. We have to use the `threshold` argument to
# only plot the simplices with diameter lower than or equal to the birth of the
# interval.

heatmap(blackhole_grayscale, aspect_ratio=1, size=(800, 800))
plot!(finite_interval, blackhole_grayscale,
      color=2,
      markershape=:d,
      markerstrokecolor=2,
      threshold=birth(finite_interval),
      label="local minimum")

plot!(infinite_interval, blackhole_grayscale,
      color=3,
      markershape=:d,
      markerstrokecolor=3,
      threshold=birth(infinite_interval),
      label="global minimum")

# Note that there are several points for each minimum because several pixels have the same
# value.

# To plot the region of the local minimum, we skip the threshold and choose a small marker
# size. We will not plot the global minimum region because it encompasses the whole image.

heatmap(blackhole_grayscale, aspect_ratio=1, size=(800, 800))
plot!(finite_interval, blackhole_grayscale,
      markersize=1,
      color=1,
      markerstrokecolor=1,
      label="local minimum region")

# We can also plot the representatives of ``H_1``, but keep in mind we have computed
# persistent cohomology. The representatives give us a chains that kill the first homology
# groups.

heatmap(blackhole_grayscale, aspect_ratio=1, size=(800, 800))
scatter!(result_min[2][1], blackhole_grayscale,
         color=1,
         label="H₁ killer1",
         markersize=0,
         linewidth=4)
scatter!(result_min[2][2], blackhole_grayscale,
         color=3,
         label="H₁ killer2",
         markersize=0,
         linewidth=4)

# Note that the ``H_1`` generators we found represent the regions around local maxima.  If
# we want to find the hole in the bright area, we must negate the image.

# Let's repeat what we just did, but with the image negated.

result_max = ripserer(Cubical(-blackhole_grayscale), cutoff=0.1, reps=true)
plot(result_max)

heatmap(blackhole_grayscale, aspect_ratio=1, size=(800, 800))
scatter!(result_max[2][1], blackhole_grayscale,
         color=1,
         label="H₁ killer",
         markersize=0,
         linewidth=4)
plot!(result_max[1][1], blackhole_grayscale,
      markershape=:d,
      color=2,
      markerstrokecolor=2,
      threshold=birth(result_max[1][1]),
      label="global maximum")
plot!(result_max[1][2], blackhole_grayscale,
      markershape=:d,
      color=3,
      markerstrokecolor=3,
      threshold=birth(result_max[1][2]),
      label="local maximum")
