# Ripserer.jl

_Flexible and efficient persistent homology computation._

Author: Matija Čufar ([@mtsch](https://github.com/mtsch/))

![](assets/title_plot.svg)

## Introduction

Ripserer is a pure Julia library for computing persistent homology based on the
[Ripser](https://github.com/Ripser/ripser) algorithm. Roughly speaking, persistent homology
detects the global topological and local geometric structure of data in a noise-resistant,
stable way. If you are unfamiliar with persistent homology, I recommend reading this
[excellent
introduction](https://towardsdatascience.com/persistent-homology-with-examples-1974d4b9c3d0).
See the Examples for further info and usage.

The main goal of this project is to provide an easy to use, generic and fast implementation
of persistent homology to the Julia ecosystem.

While this package is fully functional, it is still in development and should not be
considered stable. I try to disrupt the public interface as little as possible, but breaking
changes might still occur from time to time.

## Installation

This package is registered. To install it, simply run the following and everything should
just work.

```julia
julia> import Pkg
julia> Pkg.add("Ripserer")
```

All versions of Julia from 1.0 onward are supported, but I recommend using the latest
version of Julia for optimal performance.

## Features

Ripserer and its companion package
[PersistenceDiagrams.jl](https://github.com/mtsch/PersistenceDiagrams.jl) currently support

* Fast Vietoris-Rips and cubical persistent homology computation.
* Representative cocycle and critical simplex computation.
* Convenient persistence diagram and representative cocycle visualization via
  [Plots.jl](https://github.com/JuliaPlots/Plots.jl). [Experimental Makie.jl
  support](https://github.com/mtsch/MakieRipserer.jl) is also available.
* Bottleneck and Wasserstein matching and distance computation.
* Various persistence diagram vectorization functions, implemented with persistence images
  and persistence curves.
* Easy extensibility through a documented API.
* Computing persistent homology and representative cycles (experimental).

To access some of the features, you need to use PersistenceDiagrams.

## Performance

Much like Ripser, Ripserer uses implicit simplicial complex and reduction matrix
representations combined with the clearing optimization and other computational tricks to
achieve its speed. For a more detailed overview of these optimizations, check out [Ulrich
Bauer's article on Ripser](https://arxiv.org/abs/1908.02518).

In general, the performance of Ripserer is very close to
[Ripser](https://github.com/Ripser/ripser), usually within around 30%. Cubical homology is
up to 3× slower than that of [Cubical Ripser](https://github.com/CubicalRipser/), which uses
a more specialized algorithm. Ripserer is still a good choice for small 3d images and large
2d images. Ripserer's strength performance-wise is very sparse inputs. It also computes
some things Ripser skips, like the critical simplices. See the [Benchmarks](@ref) section for
detailed benchmarks.

## Extensibility

Ripserer is designed to be easily extended with new simplex or filtration types. See the
[Filtration Interface](@ref) and [Simplex Interface](@ref) API sections for more info. To
see an example of an extension, check out the implementation of cubical homology in
[`src/cubical.jl`](https://github.com/mtsch/Ripserer.jl/blob/master/src/cubical.j). Keep in
mind that this extension customizes almost every aspect of the algorithm and that extensions
can usually be much simpler.

If you have written an extension or are having trouble implementing one, please feel free to
open a pull request or an issue, or contact me directly.

## Contributing

All contributions are welcome, even small things like typo fixes and ideas! See the
[contribution guidelines](https://github.com/mtsch/Ripserer.jl/blob/master/CONTRIBUTING.md)
for more information.

If you used this software in a cool project, or if you have any comments, questions, or
suggestions, feel free to contact me at
[matijacufar@gmail.com](mailto:matijacufar@gmail.com).
