"""
# Ripserer.jl

Efficient computation of persistent homology.

See https://mtsch.github.io/Ripserer.jl/dev/ for documentation.
"""
module Ripserer

using Compat

using Base: @propagate_inbounds, @pure
using LinearAlgebra
using SparseArrays

using DataStructures
using Distances
using IterTools
using ProgressMeter
using RecipesBase
using TupleTools

include("infinity.jl")
export Infinity, ∞
include("primefield.jl")
export Mod

include("abstractsimplex.jl")
export AbstractSimplex, diam, coface_type, vertices, coboundary, dim
include("simplex.jl")
export IndexedSimplex, index, Simplex

include("abstractfiltration.jl")
export AbstractFiltration, n_vertices, edges, birth, threshold
include("ripsfiltration.jl")
export AbstractFlagFiltration, Rips, SparseRips

include("cubical.jl")
export Cubelet, Cubical

include("diagram.jl")
export PersistenceInterval, birth, death, persistence, representative,
    PersistenceDiagram, threshold

include("chainelement.jl")
export simplex, coefficient
include("reduction_structs.jl")
include("reduction.jl")
export ripserer
include("plotting.jl")
export barcode, barcode!


end
