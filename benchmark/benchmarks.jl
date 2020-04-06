using Ripserer
using BenchmarkTools
using SparseArrays

const SUITE = BenchmarkGroup()

include(joinpath(@__DIR__, "../test/data.jl"))

t1024 = torus(1024)
t128 = torus(128)
t64 = torus(64)

SUITE["torus 1024, dim_max = 1, modulus = 2"] =
    @benchmarkable ripserer($t1024, 1, 2)
SUITE["torus 1024, dim_max = 1, modulus = 7"] =
    @benchmarkable ripserer($t1024, 1, 7)
SUITE["torus 128,  dim_max = 2, modulus = 2"] =
    @benchmarkable ripserer($t128, 2, 2)
SUITE["torus 128,  dim_max = 2, modulus = 7"] =
    @benchmarkable ripserer($t128, 2, 7)
SUITE["torus 64,   dim_max = 4, modulus = 2"] =
    @benchmarkable ripserer($t64, 4, 2)
SUITE["torus 64,   dim_max = 4, modulus = 7"] =
    @benchmarkable ripserer($t64, 4, 7)
