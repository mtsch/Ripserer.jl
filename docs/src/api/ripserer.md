# Public API

## Ripserer

```@docs
ripserer
```

## Filtrations

```@docs
Rips
```

```@docs
SparseRips
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

## Simplices

```@docs
Simplex
```

```@docs
Cube
```

```@docs
Ripserer.vertices(::Ripserer.AbstractSimplex)
```

```@docs
Ripserer.birth(::Ripserer.AbstractSimplex)
```

```@docs
Ripserer.dim(::Ripserer.AbstractSimplex)
```

## Fields

```@docs
Mod
```

## Persistence Diagrams

See also: [PersistenceDiagrams.jl
API](https://mtsch.github.io/PersistenceDiagrams.jl/dev/api/). For convenience, Ripserer
reexports the following:

```
Ripserer.birth(::PersistenceInterval)
```

```
Ripserer.death(::PersistenceInterval)
```

```
Ripserer.persistence(::PersistenceInterval)
```

```
Ripserer.representative(::PersistenceInterval)
```

```
Ripserer.birth_simplex(::PersistenceInterval)
```

```
Ripserer.death_simplex(::PersistenceInterval)
```

```
Ripserer.barcode
```
