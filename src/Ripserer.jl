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
export AbstractFiltration, n_vertices, edges, birth
include("ripsfiltration.jl")
export AbstractFlagFiltration, RipsFiltration, SparseRipsFiltration

include("diagram.jl")
export PersistenceInterval, birth, death, persistence, representative,
    PersistenceDiagram

include("chainelement.jl")
include("reduction_structs.jl")
include("reduction.jl")
export ripserer
include("plotting.jl")
export barcode, barcode!
#=
include("cubical.jl")
=#


end
