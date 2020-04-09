module Ripserer

using Compat
using LinearAlgebra
using SparseArrays

using DataStructures
using TimerOutputs

const timer = TimerOutput()
macro timed(expr)
    esc(:(@timeit_debug timer $expr))
end
macro timed(label, expr)
    esc(:(@timeit_debug timer $label $expr))
end

include("abstract.jl")
include("rips.jl")
include("reduction.jl")

export ripserer, RipsComplex, Simplex, diam, index, coef, vertices
end
