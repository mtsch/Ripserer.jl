# # Representative Cocycles

# In this example, we will demonstrate how to compute and visualize representative cocycles
# of persistent homology generators.

# As always, we start by importing the relevant packages.
using Ripserer
using Plots
gr(); nothing # hide

# We define a function that generates some data sampled from a curve.

function curve(n)
    [(sin(t)+t/10, cos(t)+t/10) for t in range(0, 2π, length=n)]
end

plot(curve(100), legend=false, title="Curve", aspect_ratio=1, xlab="x", ylab="y")

# Then, we compute the persistent homology of that data with various numbers of points. We
# plot the representative cocycles of each result and turn them into an animation. This data
# set will always have at most one one-dimensional class. We invoke `plot(interval, data)`
# to plot the representative cocycle of the class. The same plot could be created by
# invoking `representative(interval)`. We use the `threshold_strict` argument to only plot
# the simplices with diameter strictly lower than the death time of the interval.

anim = @animate for i in vcat(3:100, 100:-1:3)
    points = curve(i)
    res = ripserer(points, reps=true)

    plt1 = plot(title="1-dimensional Representative Cocycle",
                xlab="x",
                ylab="y",
                legend=false,
                ylim=(-0.8, 1.8),
                aspect_ratio=1)
    if length(res[2]) > 0
        interval = res[2][1]
        plot!(plt1, interval, points,
              alpha=0.2,
              linewidth=2,
              threshold_strict=death(interval),
              color=3)
    end
    scatter!(plt1, points, color=2)
    plot(plt1, barcode(res), plot(res),
         layout=@layout([a [b; c]]),
         size=(800, 600))
end
mp4(anim, "cocycles_anim.mp4") # hide

# We can save the animation by running the following.

# ```julia
# mp4(anim, "cocycles_anim.mp4")
# ```
