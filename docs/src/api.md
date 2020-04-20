# API

```@docs
ripserer
```

```@docs
RipsFiltration
```

```@docs
SparseRipsFiltration
```

## Adding New Simplex Types

```@docs
AbstractSimplex
```

```@docs
index(::AbstractSimplex)
```

```@docs
coef
```

```@docs
set_coef
```

```@docs
diam(::AbstractSimplex)
```

## Adding New Filtration Types

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
SparseArrays.issparse(::AbstractFiltration)
```
