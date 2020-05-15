# # Sublevel Time Series Filtrations

# This example is based on an example from
# [ripser.py](http://ripser.scikit-tda.org/notebooks/Lower%20Star%20Time%20Series.html).

# In this example, we will show how we can use the 0-dimensional persistent homology to find
# local minima or maxima of time series data in a noise-resistant way.

# We start by loading the required packages.

using Ripserer
using Plots

using Random; gr(); Random.seed!(1337); nothing # hide

# Let's generate some data. We sample 1000 points from a noisy cosine on a slope.

n = 1000
x = range(0, 5, length=n)
y = rand(n) .+ cos.(2π * x) .+ x
plot(x, y, xlab="x", ylab="y", legend=false, title="Data")

# Even though this time series is noisy, we would like to locate the five local minima that
# can clearly be seen from the plot.

# We do this by computing 0-dimensional persistent homology of the `Cubical`. The
# resulting persistence diagram contains an interval for each local minimum. The interval's
# birth time is equal to the `y`-value at the minimum. It's death time is equal to the
# height of an adjacent local maximum. An interval with infinite persistence represents the
# global minimum.

res = ripserer(Cubical(y), dim_max=0, representatives=true)[1]
plot(res)

# We notice there is a lot of noise on the persistence diagram. We can filter it out by
# only keeping the intervals with persistence larger than 1.

res = filter(x -> persistence(x) > 1, res)
plot(res)

# After filtering we are left with the five intervals corresponding to the local minima we
# are interested in.

# We separate the infinite interval from the finite ones.

infinite = only(filter(!isfinite, res))
finite = filter(isfinite, res)
nothing # hide

# We set up a plot showing the representative of the infinite interval, which is equal to
# the whole time series. Then, we overlay a series for each of the other intervals.

plt = plot(infinite, x, y,
           seriestype=:path,
           label=string(infinite),
           legend=:bottomright,
           xlab="x", ylab="y",
           title="Persistent Local Minima")
for int in finite
    plot!(plt, int, x, y,
          seriestype=:path,
          label=string(int))
end

# Then, we add a series showing the local minima. As mentioned before, the values of the
# minima correspond to the birth times of intervals. We find the indices of the minima by
# finding the simplices with the lowest diameter in the lists of representative
# cocycles. Indexing into `x` with those will give us `x`-positions of the minima.

x_mins = x[first.(vertices.(first(sort(rep, by=diam)) for rep in representative.(res)))]
y_mins = birth.(res)
scatter!(plt, x_mins, y_mins, color=1:5, markershape=:star, label="minima")

# And finally, we plot a close-up of one of the minima to ensure we got the correct answer.

plot(res[4], x, y,
     color=4,
     seriestype=:path,
     markershape=:d,
     markersize=2,
     xlab="x", ylab="y",
     label=string(res[4]),
     title="Closeup of the 4th Local Minimum")
scatter!([x_mins[4]], [y_mins[4]], color=4, markershape=:star, label="minimum")
