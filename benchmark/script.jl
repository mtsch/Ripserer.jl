using BenchmarkTools
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Ripserer
include(joinpath(@__DIR__, "../test/data.jl"))

t1024 = torus(1024)
t256 = torus(256)
t128 = torus(128)
t64 = torus(64)

# opcija iterate_unique
bt1024_1 = @benchmark ripserer($t1024)
#=
BenchmarkTools.Trial:
  memory estimate:  497.83 MiB
  allocs estimate:  580318
  --------------
  minimum time:     40.764 s (0.28% GC)
  median time:      40.764 s (0.28% GC)
  mean time:        40.764 s (0.28% GC)
  maximum time:     40.764 s (0.28% GC)
  --------------
  samples:          1
  evals/sample:     1
=#
bt256_2 = @benchmark ripserer($t256, 2)
#=
BenchmarkTools.Trial:
  memory estimate:  1.03 GiB
  allocs estimate:  2809821
  --------------
  minimum time:     13.928 s (1.72% GC)
  median time:      13.928 s (1.72% GC)
  mean time:        13.928 s (1.72% GC)
  maximum time:     13.928 s (1.72% GC)
  --------------
  samples:          1
  evals/sample:     1
=#
bt128_3 = @benchmark ripserer($t128, 3)
#=
BenchmarkTools.Trial:
  memory estimate:  3.76 GiB
  allocs estimate:  8874208
  --------------
  minimum time:     44.526 s (1.06% GC)
  median time:      44.526 s (1.06% GC)
  mean time:        44.526 s (1.06% GC)
  maximum time:     44.526 s (1.06% GC)
  --------------
  samples:          1
  evals/sample:     1
=#
bt64_4 = @benchmark ripserer($t64, 4)
#=
BenchmarkTools.Trial:
  memory estimate:  3.96 GiB
  allocs estimate:  8370924
  --------------
  minimum time:     47.212 s (1.57% GC)
  median time:      47.212 s (1.57% GC)
  mean time:        47.212 s (1.57% GC)
  maximum time:     47.212 s (1.57% GC)
  --------------
  samples:          1
  evals/sample:     1
=#
