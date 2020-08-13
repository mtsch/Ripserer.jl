# Filtration Interface

```@docs
Ripserer.AbstractFiltration
```

```@docs
Ripserer.n_vertices(::Ripserer.AbstractFiltration)
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

# Simplex Interface

```@docs
Ripserer.AbstractSimplex
```

```@docs
dim(::Ripserer.AbstractSimplex)
```

```@docs
birth(::Ripserer.AbstractSimplex)
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
vertices(::Ripserer.AbstractSimplex)
```

```@docs
Ripserer.coboundary(::Any, ::Ripserer.AbstractSimplex)
```

```@docs
Ripserer.boundary(::Any, ::Ripserer.AbstractSimplex)
```
