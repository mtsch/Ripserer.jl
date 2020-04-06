using Ripserer
using BenchmarkTools
using SparseArrays

const SUITE = BenchmarkGroup()

include(joinpath(@__DIR__, "../test/data.jl"))

t1024 = torus(1024)
t512 = torus(512)
t256 = torus(256)
t128 = torus(128)

SUITE["torus 1024, dim_max = 1, modulus = 2"] = @benchmarkable ripserer($t1024, 1, 2), evals=1
SUITE["torus 1024, dim_max = 1, modulus = 7"] = @benchmarkable ripserer($t1024, 1, 7), evals=1
#SUITE["torus 512,  dim_max = 2, modulus = 2"] = @benchmarkable ripserer($t512, 2, 2)
#SUITE["torus 512,  dim_max = 2, modulus = 7"] = @benchmarkable ripserer($t512, 2, 7)
SUITE["torus 256,  dim_max = 3, modulus = 2"] = @benchmarkable ripserer($t256, 3, 2)
SUITE["torus 256,  dim_max = 3, modulus = 7"] = @benchmarkable ripserer($t256, 3, 7)
SUITE["torus 128,  dim_max = 4, modulus = 2"] = @benchmarkable ripserer($t128, 4, 2)
SUITE["torus 128,  dim_max = 4, modulus = 7"] = @benchmarkable ripserer($t128, 4, 7)
