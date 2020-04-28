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

```@example quick
n = 10
r = 1
R = 4
torus = [((R + r*cos(θ))*cos(φ), (R + r*cos(θ))*sin(φ), r*sin(θ))
         for θ in range(0, 2π, length=n+1)[1:end-1]
         for φ in range(0, 2π, length=n+1)[1:end-1]]
length(torus)
```

Run Ripserer.

```@example quick
using Ripserer
result = ripserer(torus)
```

Plot the result as a persistence diagram or barcode.

```@example quick
using Plots; gr()
plot(result)
savefig("diagram1.svg"); nothing # hide
barcode(result)
savefig("barcode1.svg"); nothing # hide
```

![](diagram1.svg)
![](barcode1.svg)

We notice some noise around the diagonal. This can be mitigated by running Ripserer with
`ripserer(torus, ratio=2)` or by simply filtering the diagram.

```@example quick
result[2] = filter(x -> death(x) > 2birth(x), result[2])
plot(result)
savefig("diagram2.svg"); nothing # hide
barcode(result)
savefig("barcode2.svg"); nothing # hide
```

![](diagram2.svg)
![](barcode2.svg)
