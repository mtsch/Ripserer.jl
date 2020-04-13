# simplex arithmetic microbenchmarks
module BenchArithmetic
using Ripserer
using BenchmarkTools
include(joinpath(@__DIR__, "../test/data.jl"))
suite = BenchmarkGroup()

function arithmetic(sx, α)
    res = sx
    for i in 1:1000
        λ = -coef(sx / α)
        res += sx * λ
    end
    res
end

function add(sx)
    res = sx
    for i in 1:1000
        res += sx
    end
    res
end

suite["add mod 2"] = @benchmarkable add(Simplex{2}(1,1,1))
suite["add mod 3"] = @benchmarkable add(Simplex{3}(1,1,1))
suite["add mod 7"] = @benchmarkable add(Simplex{7}(1,1,6))
suite["add mod 23"] = @benchmarkable add(Simplex{23}(1,1,6))

suite["ops mod 2"] = @benchmarkable arithmetic(Simplex{2}(2, 2, 1), 1)
suite["ops mod 3"] = @benchmarkable arithmetic(Simplex{3}(2, 2, 1), 2)
suite["ops mod 7"] = @benchmarkable arithmetic(Simplex{7}(2, 2, 3), 5)
suite["ops mod 23"] = @benchmarkable arithmetic(Simplex{23}(2, 2, 3), 5)

end

BenchArithmetic.suite
