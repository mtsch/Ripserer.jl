# Ripserer.jl

_Efficient computation of persistent homology._

![](assets/title_plot.svg)

## Introduction

Ripserer is a pure Julia library for computing persistent homology based on the
[Ripser](https://github.com/Ripser/ripser) algorithm. Roughly speaking, persistent homology
detects topological holes in data in a noise-resistant, stable way. If you are unfamiliar
with persistent homology, I recommend reading this [excellent
introduction](https://towardsdatascience.com/persistent-homology-with-examples-1974d4b9c3d0).

See the Examples for further info.

This project was created by Matija Čufar. If you used this software in your project, or if
you have any comments, questions or suggestions, feel free to contact me at
[matijacufar@gmail.com](mailto:matijacufar@gmail.com).

## Installation

This package is registered. To install it, simply run the following.

```julia
julia> using Pkg
julia> Pkg.add("Ripserer")
```

## Features

Ripserer supports the following:

* Vietoris-Rips persistent homology.
* Sublevel set persistent homology for multidimensional image and time series data.
* Calculation of persistent homology with coefficients in any (possibly user defined) field
  with the default of ``\mathbb{Z}_p`` for a prime ``p``.
* Sparse distance matrix and thresholding support.
* Computing representative cocycles of persistent cohomology.
* Plotting persistence diagrams, barcodes, matchings and representative cocycles.
* Generic API.

Ripserer uses [PersistenceDiagrams.jl](https://github.com/mtsch/PersistenceDiagrams.jl) to
represent persistence diagrams. It reexports some basic functionality, but please see that
package for more persistence diagram-related functions.

## Performance

Much like Ripser, Ripserer uses the following optimizations to achieve its speed.

* Compute persistent cohomology.
* Apply the clearing optimization.
* Don't store things that can be easily recomputed.
* Skip apparent and emergent persistence pairs.

For a detailed description of the algorithm, please see the
[Ulrich Bauer's article on Ripser](https://arxiv.org/abs/1908.02518).

In general, the performance of Ripserer is very close to
[Ripser](https://github.com/Ripser/ripser), within around 30%. Depending on the data set,
one or the other may be faster. There are no official benchmarks, because I have found
benchmarking on my computer or a CI system to be too unreliable.

## Extending

Ripserer is designed to be easily extended with new simplex or filtration types. The
interfaces are specified in the docstrings for `AbstractSimplex` and
`AbstractFiltration`. Also see the [Filtrations](@ref) and [Simplices](@ref) API sections
for more info. To see an example of an extension, check out the implementation of cubical
simplices and filtrations in
[`src/cubical.jl`](https://github.com/mtsch/Ripserer.jl/blob/master/src/cubical.j).

If you have written an extension or have trouble implementing one, please open a pull
request or an issue.

## Acknowledgments

I would like to thank:

* [@ubauer](https://github.com/ubauer) for creating the original
  [Ripser](https://github.com/Ripser/ripser) on which this project is based.
* [@ctralie](https://github.com/ctralie) and [@sauln](https://github.com/sauln) for creating
  [ripser.py](https://github.com/scikit-tda/ripser.py/) which has been a source of
  inspiration.
* Žiga Virk, for giving ideas and helping with the theoretical side of things.
