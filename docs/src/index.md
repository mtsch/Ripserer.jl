# Ripserer.jl

_Efficient computation of persistent homology._

![](assets/title_plot.svg)

## Introduction

Ripserer is a pure Julia implementation of the [Ripser](https://github.com/Ripser/ripser)
algorithm for computing persistent homology. Roughly speaking, persistent homology detects
topological holes in data in a noise-resistant, stable way. If you are unfamiliar with
persistent homology, I recommend reading this [excellent
introduction](https://towardsdatascience.com/persistent-homology-with-examples-1974d4b9c3d0).

See the Examples for further info.

## Features

Ripserer supports the following:

* Fast computation of Vietoris-Rips persistent homology.
* Calculation of persistent homology with coefficients in the field of ``\mathbb{Z}_p`` for
  a prime ``p``.
* Sparse distance matrix and thresholding support.
* Computing representative cocycles of persistent cohomology.
* Sublevel set filtrations.
* Plotting persistence diagrams and barcodes.
* Generic API.

Ripserer is not yet a full TDA software framework, since some essential features, such as
computing distances between persistence diagrams, are not available. Vietoris-Rips complexes
are also not the right fit for all kinds of data.

## Performance

Much like Ripser, Ripserer uses the following optimizations to achieve its speed.

* Compute persistent cohomology.
* Apply the clearing optimization.
* Don't store things that can be easily recomputed.
* Skip apparent and emergent persistence pairs.

For a detailed description of the algorithm, please see the
[original article](https://arxiv.org/abs/1908.02518).

In general, the performance of Ripserer is very close to
[Ripser](https://github.com/Ripser/ripser). Depending on the data set, one or the other may
be faster and the differences are usually small. There are no official benchmarks, because I
have found benchmarking on my computer or a CI system to be too unreliable.

## Extending

Ripserer is designed to be easily extended with new simplex or filtration types. There are
currently no extensions available, but implementing one should (in theory) be as simple as
overloading a few functions. The interfaces are specified in the docstrings for
`AbstractSimplex` and `AbstractFiltration`. Also see the [API](@ref) for more info.

If you have written an extension or have trouble implementing one, please open a pull
request or an issue.

## Acknowledgments

I would like to thank:

* [@ubauer](https://github.com/ubauer) for creating the original
  [Ripser](https://github.com/Ripser/ripser) on which this project is based.
* [@ctralie](https://github.com/ctralie) and [@sauln](https://github.com/sauln) for creating
  [ripser.py](https://github.com/scikit-tda/ripser.py/) which has been a source of
  inspiration.
