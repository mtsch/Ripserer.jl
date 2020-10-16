"""
# Ripserer.jl

Efficient computation of persistent homology.

See https://mtsch.github.io/Ripserer.jl/dev/ for documentation.
"""
module Ripserer

using Compat

using Base: @propagate_inbounds
using Base.Cartesian
using LinearAlgebra
using SparseArrays

using DataStructures
using Distances
using IterTools
using LightGraphs
using MiniQhull
using PersistenceDiagrams
using ProgressMeter
using RecipesBase
using StaticArrays
using TupleTools

# This functionality is imported to avoid having to deal with name clashes. There is no
# piracy involved here.
import LightGraphs: vertices, edges, nv, adjacency_matrix

import MLJModelInterface
const MMI = MLJModelInterface

# Reexporting basics makes Ripserer usable without having to import another package.
import PersistenceDiagrams: birth, threshold, dim
export birth, death, persistence, representative, birth_simplex, death_simplex, barcode

export Mod,
    Simplex,
    Cube,
    index,
    vertices,
    dim,
    threshold,
    simplex,
    coefficient,
    Rips,
    SparseRips,
    Cubical,
    Custom,
    Alpha,
    ripserer,
    reconstruct_cycle

include("primefield.jl")
include("simplex.jl")
include("abstractfiltration.jl")
include("chainelement.jl")

include("computation/utils.jl")
include("computation/zerodimensional.jl")
include("computation/workingchain.jl")
include("computation/reducedmatrix.jl")
include("computation/coboundarymatrices.jl")
include("computation/ripserer.jl")

include("filtrations/utils.jl")
include("filtrations/rips.jl")
include("filtrations/cubical.jl")
include("filtrations/custom.jl")
include("filtrations/alpha.jl")

include("cycles.jl")

include("mlj.jl")
include("simplexrecipes.jl")

end
