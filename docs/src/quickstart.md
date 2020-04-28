```@setup layouts
using Plots; gr()
```

# Quick Start

This package is still under development and is currently unregistered. To install it, run
the following.

```julia
using Pkg
Pkg.add("https://github.com/mtsch/Ripserer.jl")
```

Generate 100 points sampled from a torus.

```@example data
n = 10
r = 1
R = 4
torus = [((R + r*cos(θ))*cos(φ), (R + r*cos(θ))*sin(φ), r*sin(θ))
         for θ in range(0, 2π, length=n+1)[1:end-1]
         for φ in range(0, 2π, length=n+1)[1:end-1]]
```

Run Ripserer.

```@example run
using Ripserer
result = ripserer(torus)
```

Plot the result as a persistence diagram or barcode.

```@example plot
julia> using Plots; gr()
julia> plot(result)
```

```@example plot
julia> barcode(result)
```

We notice some noise around the diagonal. This can be mitigated by running Ripserer with
`ripserer(torus, ratio=2)` or by simply filtering the diagram.

```@example filter
result[2] = filter(x -> death(x) > 2birth(x), result[2])
```
```@example filter
julia> plot(result)
```
```@example filter
julia> barcode(result)
```
