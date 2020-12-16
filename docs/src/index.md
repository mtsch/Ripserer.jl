![](assets/logo-title.svg)

_Flexible and efficient persistent homology computation._

Author: Matija Čufar ([@mtsch](https://github.com/mtsch/))

## Introduction

Ripserer is a pure Julia library for computing persistent homology based on the
[Ripser](https://github.com/Ripser/ripser) algorithm. Roughly speaking, persistent homology
detects the global topological and local geometric structure of data in a noise-resistant,
stable way. If you are unfamiliar with persistent homology, I recommend reading this
[excellent
introduction](https://towardsdatascience.com/persistent-homology-with-examples-1974d4b9c3d0).

Please see the [Usage Guide](@ref) for a quick introduction, and the [API](@ref) page for
detailed descriptions of Ripserer's functionality.

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

* Fast Vietoris-Rips and cubical, and alpha complex persistent homology computation.
* Representative cocycle, cycle, and critical simplex computation.
* Convenient persistence diagram and representative cocycle visualization via
  [Plots.jl](https://github.com/JuliaPlots/Plots.jl). Experimental
  [Makie.jl](https://github.com/JuliaPlots/Makie.jl) is also available
  [here](https://github.com/mtsch/MakieRipserer.jl).
* Bottleneck and Wasserstein matching and distance computation.
* Various persistence diagram vectorization functions, implemented with persistence images
  and persistence curves.
* Easy extensibility through a documented API.
* Experimental shortest representative cycle computation.
* Experimental sparse circular coordinate computation.

To access some of the features, you need to use the PersistenceDiagrams.jl package.

## Performance

Much like Ripser, Ripserer uses several computational tricks to achieve its speed. Among
others, these include an implicit simplicial complex representation and the clearing
optimization. For a more detailed overview of these optimizations, check out [Ulrich Bauer's
article on Ripser](https://arxiv.org/abs/1908.02518).

In general, the performance of Ripserer is very close to
[Ripser](https://github.com/Ripser/ripser), usually within around 30%. Ripserer's strength
performance-wise is very sparse inputs, where it can sometimes outperform Ripser. It also
computes some things Ripser skips, like the critical simplices.

Ripserer's Cubical homology is up to 3× slower than that of [Cubical
Ripser](https://github.com/CubicalRipser/), which uses a more specialized
algorithm. Ripserer is still a good choice for small 3d images and large 2d images. Unlike
Cubical Ripser, it also supports computations on images of dimensions higher than 4.

See the [Benchmarks](@ref) section for more detailed benchmarks.

## Extending

Ripserer is designed to be easily extended with new simplex or filtration types. See the
[Abstract Types and Interfaces](@ref) API section for more information.

If you have written an extension or are having trouble implementing one, please feel free to
open a pull request or an issue. You may also contact me directly.

## Contributing

All contributions are welcome, even small things like typo fixes and ideas! See the
[contribution guidelines](https://github.com/mtsch/Ripserer.jl/blob/master/CONTRIBUTING.md)
for more information.

If you used this software in a cool project, or if you have any comments, questions, or
suggestions, feel free to contact me at
[matijacufar@gmail.com](mailto:matijacufar@gmail.com).

## Citing

If you used Ripserer in your work, consider citing the [JOSS
paper](https://joss.theoj.org/papers/10.21105/joss.02614).

A bibtex entry is provided in
[CITATION.bib](https://github.com/mtsch/Ripserer.jl/blob/master/CITATION.bib).
