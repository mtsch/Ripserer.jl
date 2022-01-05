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
using Graphs
using MiniQhull
using PersistenceDiagrams
using ProgressMeter
using RecipesBase
using StaticArrays
using TupleTools

import MLJModelInterface

# This functionality is imported to avoid having to deal with name clashes. There is no
# piracy involved here.
import Graphs: vertices, edges, nv, adjacency_matrix

# Reexporting basics makes Ripserer usable without having to import another package.
import PersistenceDiagrams: birth, threshold, dim
export birth,
    death, persistence, midlife, representative, birth_simplex, death_simplex, barcode

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
    Cubical,
    Custom,
    Alpha,
    EdgeCollapsedRips,
    ripserer,
    reconstruct_cycle,
    Partition,
    CircularCoordinates,
    RipsPersistentHomology,
    AlphaPersistentHomology,
    CubicalPersistentHomology

include("base/primefield.jl")
include("base/abstractcell.jl")
include("base/abstractfiltration.jl")
include("base/abstractsimplex.jl")
include("base/simplex.jl")
include("base/chainelement.jl")
include("base/chain.jl")
include("base/simplexrecipes.jl")

include("computation/utils.jl")
include("computation/zerodimensional.jl")
include("computation/reducedmatrix.jl")
include("computation/cohomology.jl")
include("computation/ripserer.jl")

include("filtrations/utils.jl")
include("filtrations/rips.jl")
include("filtrations/cubical.jl")
include("filtrations/custom.jl")
include("filtrations/alpha.jl")
include("filtrations/edgecollapse.jl")

include("extra/cycles.jl")
include("extra/circularcoordinates.jl")
include("extra/mlj.jl")

end
