# # Stability

# In this example we will demonstrate the stability of persistent homology. The stability
# theorem roughly states, that a small change in the input data will result in a small
# change in the resulting persistence diagram. In other words, persistent homology is very
# tolerant of noisy data. Also see the
# [Distances](https://mtsch.github.io/PersistenceDiagrams.jl/dev/generated/distances/)
# example in [PersistenceDiagrams.jl](https://github.com/mtsch/PersistenceDiagrams.jl).

# Again, start with loading some packages.

using Ripserer
using Plots
using Random; Random.seed!(1337); gr(); nothing # hide

# As in the Basics example, we will look at the persistent homology of a noisy circle.

function noisy_circle(n; r1=1, r2=1, noise=0.1)
    points = NTuple{2, Float64}[]
    for _ in 1:n
        θ = 2π * rand()
        push!(points, (r1*sin(θ) + noise*rand() - noise/2,
                       r2*cos(θ) + noise*rand() - noise/2))
    end
    points
end

# We will look at the first persistent homology group of this space.

# First, let's see what happens if we repeatedly sample 100 random points from a circle.

anim = @animate for _ in 1:200
    points = noisy_circle(100, noise=0)
    result = ripserer(points)

    plt_pts = scatter(points, legend=false,
                      aspect_ratio=1,
                      xlim=(-2.2, 2.2),
                      ylim=(-2.2, 2.2),
                      title="Data")
    plt_diag = plot(result, infinity=3)

    plot(plt_pts, plt_diag, size=(800, 400))
end
gif(anim, "stability_anim_1.gif") # hide

# We notice that an interval in ``H_1`` always stands out and that its death remains
# constant. The only thing that changes is the birth time. The birth time is equal to the
# largest distance between adjacent points in the circle. At birth time, the circle is
# connected.

# Now, let's add some noise!

anim = @animate for _ in 1:200
    points = noisy_circle(100, noise=0.2)
    result = ripserer(points)

    plt_pts = scatter(points, legend=false,
                      aspect_ratio=1,
                      xlim=(-2.2, 2.2),
                      ylim=(-2.2, 2.2),
                      title="Data")
    plt_diag = plot(result, infinity=3)

    plot(plt_pts, plt_diag, size=(800, 400))
end
gif(anim, "stability_anim_2.gif") # hide

# The interval is jumping around a lot more now, but it hovers around the same general area.
# It's still clearly the most persistent feature of our space.

# Next, let's look at how adding more and more noise affects the diagram.

anim = @animate for noise in vcat(0:0.01:1, 1:-0.01:0)
    points = noisy_circle(100, noise=noise)
    result = ripserer(points)

    plt_pts = scatter(points, legend=false,
                      aspect_ratio=1,
                      xlim=(-2.2, 2.2),
                      ylim=(-2.2, 2.2),
                      title="Data")
    plt_diag = plot(result, infinity=3)

    plot(plt_pts, plt_diag, size=(800, 400))
end
gif(anim, "stability_anim_3.gif") # hide

# We see we have to add quite a bit of noise to destroy the diagram. Notice how the death
# time of the interval decreases as we add noise. This is the result of the diameter of the
# hole in our circle shrinking.

# Finally, let's stretch our circle.

anim = @animate for r in vcat(0.0:0.02:2, 2:-0.02:0.0)
    points = noisy_circle(100, noise=0.1, r1=r)
    result = ripserer(points)

    plt_pts = scatter(points, legend=false,
                      aspect_ratio=1,
                      xlim=(-2.2, 2.2),
                      ylim=(-2.2, 2.2),
                      title="Data")
    plt_diag = plot(result, infinity=3)

    plot(plt_pts, plt_diag, size=(800, 400))
end
gif(anim, "stability_anim_4.gif") # hide

# Again, we see the persistent homology stays stable, as long as the data at least somewhat
# resembles a circle.
