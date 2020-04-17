module Ripserer

using Compat
using Future: copy!

using Base: @propagate_inbounds

using LinearAlgebra
using SparseArrays

using DataStructures

include("interface.jl")
include("abstract.jl")
include("simplex.jl")
include("rips.jl")
include("reduction.jl")

export ripserer,
    AbstractFiltration, RipsFiltration, SparseRipsFiltration, AbstractSimplex, Simplex,
    coef, set_coef, index, diam,
    vertices, dim_max, threshold, dist, edges
end
