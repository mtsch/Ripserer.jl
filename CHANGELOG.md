# v0.16.9

Replace LightGraps with Graphs.

# v0.16.4

Add [MLJ.jl](https://github.com/alan-turing-institute/MLJ.jl) support.

# v0.16.3

New function: `midlife`.

# v0.16.0

External changes:

* `progress` keyword argument renamed to `verbose`, `field_type` keyword argument renamed to
  `field`.
* New interface: `ripserer(::Type{AbstractFiltration}, args...; kwargs...)`.
* Added `CircularCoordinates`.

Interface changes:

* All filtration constructors now have to take `verbose` as a keyword argument.
* Replaced vectors of `ChainElement`s with `Chain`s.
* Added `AbstractCell`.
* `Cube` is now an `AbstractCell`, `AbstractSimplex` is reserved for actual simplices.
* Simplices are no longer `Array`s.
* `simplex`, `unsafe_simplex`, and `unsafe_cofacet` no longer take a `sign` argument.

# v0.15.4

* Use `PersistenceDiagrams` v0.8.

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
* Involuted homology can be computed with the `alg=:involuted` keyword argument.
* The `reps` keyword argument can be set to a collection of integers, finding
  representatives only for specified dimensions.
* Implicit or explicit reduction can be set with the `implicit` keyword argument.
* Improved progress printing.
* New function: `find_apparent_pairs`.
