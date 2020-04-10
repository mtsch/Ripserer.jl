using BenchmarkTools
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Ripserer
using Ripser
include(joinpath(@__DIR__, "../test/data.jl"))

t1024 = torus(1024)
t256 = torus(256)
t128 = torus(128)
t64 = torus(64)

# opcija iterate_unique
bt1024_1_ripserer = @benchmark ripserer($t1024)
bt1024_1_ripser = @benchmark ripser($t1024)

bt256_2_ripserer = @benchmark ripserer($t256, dim_max=2)
bt256_2_ripser = @benchmark ripser($t256, dim_max=2)

bt128_3_ripserer = @benchmark ripserer($t128, dim_max=3)
bt128_3_ripser = @benchmark ripser($t128, dim_max=3)
# 67

bt64_4_ripserer = @benchmark ripserer($t64, dim_max=4)
bt64_4_ripser = @benchmark ripser($t64, dim_max=4)
