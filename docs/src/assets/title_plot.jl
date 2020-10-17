# This file generates the plot in the landing page:w
using Ripserer
using Plots;
gr();
using Random;
Random.seed!(1337);

function noisy_circle(n; r=0.9, noise=0.1, center)
    points = NTuple{2,Float64}[]
    x, y = center
    for _ in 1:n
        θ = 2π * rand()
        push!(
            points,
            (
                x + r * sin(θ) + noise * rand() - noise / 2,
                y + r * cos(θ) + noise * rand() - noise / 2,
            ),
        )
    end
    return points
end

data = [
    noisy_circle(200; noise=0.1, center=(-1, 0))
    noisy_circle(200; noise=0.2, center=(1, 0))
    noisy_circle(200; noise=0.3, center=(0, √3))
]

res = ripserer(data; reps=true)
plt_dgm = plot(res; markeralpha=1, markersize=3)
plt_data = plot(;
    legend=false,
    aspect_ratio=1,
    title="Representative Cocycles",
    xlim=(-2.1, 2.1),
    ylim=(-1.1, 3.1),
    xlabel=:x,
    ylabel=:y,
)

for (i, int) in enumerate(sort(res[2]; by=death, rev=true))
    if i == 1
        col = 2 # red
    elseif i == 2
        col = 4 # purple
    elseif i == 3
        col = 3 # green
    else
        col = :black
    end

    plot!(plt_data, int, data; alpha=0.1, color=col)
end
scatter!(plt_data, data; markersize=1)
plt_data

plot(plt_data, plt_dgm; size=(800, 400))
savefig("title_plot.svg")
