# API

## Ripserer

```@docs
ripserer
```

```@docs
RipsFiltration
```

```@docs
SparseRipsFiltration
```

## Persistence Intervals and Diagrams

```@docs
PersistenceInterval
```

```@docs
PersistenceDiagram
```

```@docs
RecipesBase.plot(::Union{PersistenceDiagram, AbstractVector{<:PersistenceDiagram}})
```

```@docs
barcode(::Union{PersistenceDiagram, AbstractVector{<:PersistenceDiagram}})
```

## Simplex Types

```@docs
AbstractSimplex
```

```@docs
index(::AbstractSimplex)
```

```@docs
coef(::AbstractSimplex)
```

```@docs
set_coef(::AbstractSimplex, ::Any)
```

```@docs
diam(::AbstractSimplex)
```

```@docs
dim(::AbstractSimplex)
```

```@docs
vertices(::AbstractSimplex)
```

```@docs
coface_type(::AbstractSimplex)
```

```@docs
coboundary
```

## Filtration Types

```@docs
AbstractFiltration
```

```@docs
AbstractFlagFiltration
```

```@docs
n_vertices
```

```@docs
edges
```

```@docs
diam(::AbstractFiltration, ::Any)
```

```@docs
diam(::AbstractFiltration, ::Any, ::Any, ::Any)
```

```@docs
SparseArrays.issparse(::AbstractFiltration)
```
