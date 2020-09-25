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

## Experimental Features

```@docs
reconstruct_cycle
```
