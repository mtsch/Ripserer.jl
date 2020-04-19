module Ripserer

using Compat
using Future: copy!

using Base: @propagate_inbounds
using LinearAlgebra
using SparseArrays

using DataStructures
using Distances

include("convenience.jl")
include("interface.jl")
include("simplex.jl")
include("rips.jl")
include("coboundary.jl")
include("reduction.jl")

export AbstractFiltration, AbstractSimplex,
    Simplex, coef, set_coef, index, diam,
    AbstractFlagFiltration, RipsFiltration, SparseRipsFiltration,
    n_vertices, threshold, edges,
    ripserer

end
