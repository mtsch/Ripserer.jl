# Ripserer.jl

[![Coverage Status](https://coveralls.io/repos/github/mtsch/Ripserer.jl/badge.svg?branch=master)](https://coveralls.io/github/mtsch/Ripserer.jl?branch=master)
[![codecov](https://codecov.io/gh/mtsch/Ripserer.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mtsch/Ripserer.jl)
[![Build Status](https://travis-ci.org/mtsch/Ripserer.jl.svg?branch=master)](https://travis-ci.org/mtsch/Ripserer.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/cc709npw3lp76yc8?svg=true)](https://ci.appveyor.com/project/mtsch/ripserer-jl)
[![Aqua QA](https://img.shields.io/badge/Aqua.jl-%F0%9F%8C%A2-aqua.svg)](https://github.com/tkf/Aqua.jl)

A Julia reimplementation of the [ripser](https://github.com/Ripser/ripser) algorithm for
persistent homology. This package is not a direct translation and might do or name some
things differently.

Ripserer's performance is generally around 2 times slower than
[ripser](https://github.com/Ripser/ripser), but in some cases, it performs just as well or
even better.

The goal of this project is to learn how [ripser](https://github.com/Ripser/ripser) works
and to make it easy to play with extensions.

## Usage

```julia
    ripserer(dists::AbstractMatrix{T}; dim_max=1, modulus=2, threshold)
    ripserer(filtration::AbstractFiltration)
```

Compute the persistent homology of metric space represented by distance matrix `dists` or
`filtration`. If `dist` is sparse, a sparse filtration is constructed.

Keyword Arguments:

* `dim_max`: compute persistent homology up to this dimension.
* `modulus`: compute persistent homology with coefficients in the prime field of integers
             mod `modulus`.
* `threshold`: compute persistent homology up to diameter smaller than threshold.
               Defaults to radius of input space.

Returns a vector of persistence diagrams, each represented by a vector of `(birth, death)`
where `birth` and `death` are of the type `eltype(dists)`.

## Extending

One of the goals of this project is to make it easy to extend with new filtration and
simplex types. A subtype of `AbstractFiltration` or `AbstractSimplex` should only _need_ to
implement the interface functions found in [`src/interface.jl`][src/interface.jl]. Other
functions may need to be overloaded to improve performance.
