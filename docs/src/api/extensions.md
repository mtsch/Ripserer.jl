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
Ripserer.birth(::Ripserer.AbstractFiltration, ::Any)
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
Ripserer.AbstractRipsFiltration
```

```@docs
Ripserer.dist
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
index(::Ripserer.AbstractSimplex)
```

```@docs
Base.sign(::Ripserer.AbstractSimplex)
```

```@docs
Base.:-(::Ripserer.AbstractSimplex)
```

```@docs
Ripserer.coboundary(::Any, ::Ripserer.AbstractSimplex)
```

```@docs
Ripserer.boundary(::Any, ::Ripserer.AbstractSimplex)
```
