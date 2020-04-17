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
length(::AbstractFiltration)
```

```@docs
dist
```

```@docs
edges
```

```@docs
dim_max
```

```@docs
diam(::AbstractFiltration, ::Any)
```

```@docs
binomial(::AbstractFiltration, ::Any, ::Any)
```

```@docs
threshold
```

```@docs
SparseArrays.issparse(::AbstractFiltration)
```
