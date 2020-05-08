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

include("simplices/abstractsimplex.jl")
export AbstractSimplex, diam, coface_type, vertices, coboundary, dim
include("simplices/indexedsimplex.jl")
export IndexedSimplex, index
include("simplices/simplex.jl")
export Simplex

include("filtrations/abstractfiltration.jl")
export AbstractFiltration, n_vertices, edges, birth
include("filtrations/ripsfiltration.jl")
export AbstractFlagFiltration, RipsFiltration, SparseRipsFiltration

include("primefield.jl")
export PrimeField

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
