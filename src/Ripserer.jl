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
using TupleTools

# reexporting basics these makes Ripserer usable without having to import another package.
import PersistenceDiagrams: birth, threshold, dim
export
    birth, death, persistence, representative, birth_simplex, death_simplex, barcode

export
    Mod, Simplex, Cubelet, index, vertices, dim, threshold, simplex, coefficient,
    Rips, SparseRips, Cubical, ripserer


include("primefield.jl")

include("simplex.jl")

include("abstractfiltration.jl")
include("ripsfiltration.jl")
include("cubical.jl")

include("chainelement.jl")
include("zerodimensional.jl")
include("reductionmatrix.jl")
include("main.jl")

include("simplexrecipes.jl")

end
