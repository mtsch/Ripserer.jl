# v0.15.3

* Fix type instability in zeroth interval generation.

# v0.15.2

* Update compat with Distances.jl.

# v0.15.1

* Representative cocycles are computed for infinite intervals.

# v0.15.0

* `SparseRips(...)` is deprecated. Use `Rips(...; sparse=true)`.
* `AbstractFiltration`s now need to define `births` instead of `birth`.
* Results are now sorted by persistence instead of birth time.
* Homology is now computed with the `alg=:homology` keyword argument.
* Assisted homology can be computed with the `alg=:assisted` keyword argument.
* The `reps` keyword argument can be set to a collection of integers, finding
  representatives only for specified dimensions.
* Implicit or explicit reduction can be set with the `implicit` keyword argument.
* Improved progress printing.
* New function: `find_apparent_pairs`.
