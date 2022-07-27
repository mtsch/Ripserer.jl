# Related Julia Packages

This section attempts to provide an overview of Julia packages implementing persistent
homology. If you're trying to do something Ripserer is missing, one of these might have what
you're looking for.

This list is incomplete and only includes things I'm aware of. The descriptions are based
on my (limited) experience with these packages.

* [ComputationalHomology.jl](https://github.com/wildart/ComputationalHomology.jl) uses a
  different, slower algorithm, but can do some things Ripserer can't such as Čech
  persistent homology and homology of CW complexes. The package is a part of
  [TDA.jl](https://github.com/wildart/TDA.jl), which also offers other topological data
  analysis tools like Mapper.

* [Eirene.jl](https://github.com/Eetion/Eirene.jl) uses a [different algorithm based on
  matroids](https://arxiv.org/abs/1606.00199). A benefit of this algorithm is that it can
  recover persistent homology generators and a bunch of other data.

* [Ripser.jl](https://github.com/mtsch/Ripser.jl) my deprecated wrapper of the original
  C++ program. Very bare-bones and outdated. Runs a bit slower than Ripserer.

* [Sparips.jl](https://github.com/bbrehm/Sparips.jl) this is a preprocessor that comes with
  a different wrapper of Ripser. The preprocessor allows you to compute persistent homology
  of very large datasets. Integrating it with Ripserer is on my TODO list.

* [PersistentCohomology.jl](https://github.com/piever/PersistentCohomology.jl) does
  essentially the same thing, but with a different algorithm, that is orders of magnitude
  slower.

If you are a developer of a persistent homology package not on this list, or any of the
information here is incorrect, let me know or open a PR.
