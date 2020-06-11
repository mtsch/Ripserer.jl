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
using StaticArrays

# reexporting basics these makes Ripserer usable without having to import another package.
import PersistenceDiagrams: birth, threshold
export PersistenceDiagram, PersistenceInterval
export birth, death, persistence, representative, barcode

export Mod
export AbstractSimplex, diam, coface_type, face_type, vertices, coboundary, boundary, dim
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
include("zerodimensional.jl")
include("reductionmatrix.jl")
include("ripserer.jl")

include("simplexrecipes.jl")

end
