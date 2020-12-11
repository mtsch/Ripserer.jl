# API

## Ripserer

```@docs
ripserer
```

## Filtrations

```@docs
Rips
```

```@docs
Cubical
```

```@docs
Custom
```

```@docs
Alpha
```

```@docs
EdgeCollapsedRips
```

## Persistence Diagrams

Persistence diagrams live in a separate package,
[PersistenceDiagrams.jl](https://github.com/mtsch/PersistenceDiagrams.jl). The package is
documented in detail [here](https://mtsch.github.io/PersistenceDiagrams.jl/dev/).

If you are looking for
Wasserstein or bottleneck distances, persistence images, betti curves, landscapes, and
similar, you will need to run `using PersistenceDiagrams`.

For convenience, the following basic functionality is reexported by Ripserer:

```@docs
PersistenceDiagrams.PersistenceDiagram
```

```@docs
PersistenceDiagrams.PersistenceInterval
```

```@docs
birth(::PersistenceDiagrams.PersistenceInterval)
```

```@docs
death(::PersistenceDiagrams.PersistenceInterval)
```

```@docs
persistence(::PersistenceDiagrams.PersistenceInterval)
```

```@docs
midlife
```

```@docs
representative(::PersistenceDiagrams.PersistenceInterval)
```

```@docs
birth_simplex(::PersistenceDiagrams.PersistenceInterval)
```

```@docs
death_simplex(::PersistenceDiagrams.PersistenceInterval)
```

```@docs
barcode
```

## Simplices and Representatives

```@docs
Ripserer.Simplex
```

```@docs
Ripserer.Cube
```

```@docs
Ripserer.dim(::Ripserer.AbstractCell)
```

```@docs
Ripserer.birth(::Ripserer.AbstractCell)
```

```@docs
Ripserer.index(::Ripserer.AbstractCell)
```

```@docs
Ripserer.vertices(::Ripserer.AbstractCell)
```

```@docs
Ripserer.Chain
```

```@docs
Mod
```

## Experimental Features

```@docs
Ripserer.reconstruct_cycle
```

```@docs
Ripserer.Partition
```

```@docs
Ripserer.CircularCoordinates
```

## Abstract Types and Interfaces

```@docs
Ripserer.AbstractFiltration
```

```@docs
Ripserer.nv(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.births(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.vertices(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.edges(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.simplex_type
```

```@docs
Ripserer.simplex
```

```@docs
Ripserer.unsafe_simplex
```

```@docs
Ripserer.unsafe_cofacet
```

```@docs
Ripserer.threshold(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.columns_to_reduce
```

```@docs
Ripserer.emergent_pairs
```

```@docs
Ripserer.postprocess_diagram
```

```@docs
Ripserer.distance_matrix
```

```@docs
Ripserer.AbstractRipsFiltration
```

```@docs
Ripserer.adjacency_matrix(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.AbstractCustomFiltration
```

```@docs
Ripserer.simplex_dicts
```

```@docs
Ripserer.AbstractCell
```

```@docs
Ripserer.AbstractSimplex
```

```@docs
Base.sign(::Ripserer.AbstractCell)
```

```@docs
Base.:-(::Ripserer.AbstractCell)
```

```@docs
Ripserer.coboundary
```

```@docs
Ripserer.boundary
```
