# Interfaces

## Filtrations

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

## Simplices

```@docs
Ripserer.AbstractCell
```

```@docs
Ripserer.AbstractSimplex
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

```@docs
Ripserer.Simplex
```

```@docs
Ripserer.Cube
```
