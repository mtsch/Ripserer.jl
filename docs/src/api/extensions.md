# Interfaces

## Filtration Interface

```@docs
Ripserer.AbstractFiltration
```

```@docs
Ripserer.nv(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.edges(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.births(::Ripserer.AbstractFiltration)
```

```@docs
Ripserer.threshold(::Ripserer.AbstractFiltration)
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
Ripserer.AbstractCustomFiltration
```

```@docs
Ripserer.simplex_dicts
```

```@docs
Ripserer.adjacency_matrix(::Ripserer.AbstractFiltration)
```

## Simplex Interface

```@docs
Ripserer.AbstractSimplex
```

```@docs
index(::Ripserer.AbstractCell)
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
