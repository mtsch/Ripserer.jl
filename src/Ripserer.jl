module Ripserer

using Compat

using Base: @propagate_inbounds, @pure
using LinearAlgebra
using SparseArrays

using DataStructures
using Distances
using RecipesBase
using TupleTools

include("infinity.jl")
export Infinity, ∞
include("abstractsimplex.jl")
export AbstractSimplex, coef, set_coef, index, diam, coface_type, vertices, coboundary, dim
include("simplex.jl")
export Simplex
include("abstractfiltration.jl")
export AbstractFiltration, n_vertices, edges, birth, issparse
include("ripsfiltration.jl")
export AbstractFlagFiltration, RipsFiltration, SparseRipsFiltration
include("diagram.jl")
export PersistenceInterval, birth, death, persistence, representative,
    PersistenceDiagram, barcode, barcode!
include("reduction.jl")
export ripserer

end
