module Ripserer

using Compat
using Future: copy!

using Base: @propagate_inbounds, @pure
using LinearAlgebra
using SparseArrays

using DataStructures
using Distances
using RecipesBase
using TupleTools

include("interface.jl")
include("simplex.jl")
include("rips.jl")
include("diagrams.jl")
include("reduction.jl")

export Infinity, ∞,
    AbstractFiltration, AbstractSimplex,
    Simplex, coef, set_coef, index, diam, coface_type, vertices, coboundary, dim,
    AbstractFlagFiltration, RipsFiltration, SparseRipsFiltration,
    n_vertices, threshold, edges, edge_type,
    PersistenceInterval, birth, death, cocycle, PersistenceDiagram,
    barcode, barcode!,
    ripserer

end
