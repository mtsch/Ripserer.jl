module Ripserer

using Compat
using LinearAlgebra
using SparseArrays

using DataStructures

include("abstract.jl")
include("rips.jl")
include("reduction.jl")

export ripserer,
    AbstractFiltration, RipsFiltration, SparseRipsFiltration, AbstractSimplex, Simplex,
    diam, index, coef, vertices, dim_max, threshold, dist
end
