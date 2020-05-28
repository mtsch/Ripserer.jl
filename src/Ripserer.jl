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
using Hungarian
using IterTools
using PersistenceDiagrams
using ProgressMeter
using RecipesBase
using TupleTools

# reexport basics from PersistenceDiagrams
import PersistenceDiagrams: birth
export PersistenceDiagram, PersistenceInterval
export birth, death, persistence, representative, barcode, barcode!

export Mod
export AbstractSimplex, diam, coface_type, vertices, coboundary, dim
export IndexedSimplex, index, Simplex
export AbstractFiltration, n_vertices, edges, threshold
export AbstractFlagFiltration, Rips, SparseRips
export Cubelet, Cubical
export simplex, coefficient
export ripserer

include("primefield.jl")

include("abstractsimplex.jl")
include("simplex.jl")

include("abstractfiltration.jl")
include("ripsfiltration.jl")
include("cubical.jl")

include("chainelement.jl")
include("reduction_structs.jl")
include("reduction.jl")
include("simplexrecipes.jl")

end
