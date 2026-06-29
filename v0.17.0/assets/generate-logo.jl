# This file generates the plot in the logo
using Ripserer
using Plots
gr()
using Random
Random.seed!(1337)

function circle(n; r=0.7, center, s=0)
    points = NTuple{2,Float64}[]
    x, y = center
    for θ in range(0, 2π; length=n + 1)[circshift(1:n, s)]
        push!(points, (x + r * sin(θ), y + r * cos(θ)))
    end
    return points
end

data = [
    circle(20; center=(-1, 0), s=0)
    circle(20; center=(1, 0), s=0)
    circle(20; center=(0, √3), s=0)
]

res = ripserer(data; reps=true, dim_max=2)

ripserer_logo = plot(;
    aspect_ratio=1,
    xlim=(-2.1, 2.1),
    ylim=(-1.1, 3.1),
    xlabel=:x,
    ylabel=:y,
    ticks=nothing,
    legend=false,
    border=:none,
)
for (i, int) in enumerate(sort(res[2]; by=death, rev=true))
    if i == 1
        col = 2 # red
        α = 1
    elseif i == 2
        col = 4 # purple
        α = 1
    elseif i == 3
        col = 3 # green
        α = 1
    else
        col = 1
        α = 0.2
    end

    plot!(
        ripserer_logo,
        int,
        data;
        alpha=α,
        color=col,
        threshold_strict=death(int),
        linewidth=5,
    )
end
scatter!(ripserer_logo, data; markersize=4, color=1)
display(ripserer_logo)

diagrams_logo = plot(
    res;
    title="",
    ticks=nothing,
    legend=false,
    xlabel="",
    ylabel="",
    markersize=3,
    size=(100, 100),
)
