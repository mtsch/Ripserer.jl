# Ripserer.jl

_Efficient computation of persistent homology._

A Julia implementation of the [ripser](https://github.com/Ripser/ripser) algorithm for
persistent homology.

Like [ripser](https://github.com/Ripser/ripser), Ripserer uses the following optimizations
to achieve fast computation of Vietoris-Rips persistent homology.

* Compute persistent cohomology
* Apply the clearing optimization
* Don't store things that can be easily recomputed
* Skip apparent and emergent persistence pairs.

For a detailed description of the algorithm, please see the
[original article](https://arxiv.org/abs/1908.02518).

## Performance

In general, the performance of Ripserer is very close to
[ripser](https://github.com/Ripser/ripser). Depending on the data set, one or the other may
be faster and the differences are usually small.

## Extending

Ripserer is designed to be easily extended with new simplex or filtration types. There are
currently no extensions available, but implementing one should (in theory) be as simple as
overloading a few functions. The interfaces are specified in the docstrings for
`AbstractSimplex` and `AbstractFiltration`.

## Manual

```@contents
Pages = ["quickstart.md",
         "api.md",
        ]
Depth = 1
```
