# API

## Ripserer

```@docs
ripserer
```

## Persistence Intervals and Diagrams

```@docs
PersistenceInterval
```

```@docs
birth(::PersistenceInterval)
```

```@docs
death(::PersistenceInterval)
```

```@docs
persistence(::PersistenceInterval)
```

```@docs
representative(::PersistenceInterval)
```

```@docs
PersistenceDiagram
```

```@docs
barcode(::Union{PersistenceDiagram, AbstractVector{<:PersistenceDiagram}})
```

## Filtration Types

```@docs
AbstractFiltration
```

```@docs
n_vertices(::AbstractFiltration)
```

```@docs
edges(::AbstractFiltration)
```

```@docs
diam(::AbstractFiltration, ::Any)
```

```@docs
diam(::AbstractFiltration, ::AbstractSimplex, ::Any, ::Any)
```

```@docs
birth(::AbstractFiltration, ::Any)
```

```@docs
threshold(::AbstractFiltration)
```

```@docs
AbstractFlagFiltration
```

```@docs
Rips
```

```@docs
SparseRips
```

```@docs
Cubical
```

## Simplex Types

```@docs
AbstractSimplex
```

```@docs
dim(::AbstractSimplex)
```

```@docs
diam(::AbstractSimplex)
```

```@docs
Base.sign(::AbstractSimplex)
```

```@docs
Base.:-(::AbstractSimplex)
```

```@docs
coface_type(::AbstractSimplex)
```

```@docs
vertices(::AbstractSimplex)
```

```@docs
coboundary(::Any, ::AbstractSimplex)
```

```@docs
IndexedSimplex
```

```@docs
index(::AbstractSimplex)
```

```@docs
Simplex
```

```@docs
Cubelet
```
