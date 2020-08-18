using Ripserer
using Test

@testset "Custom filtration" begin
    flt = Custom([
        (1,) => 0,
        (4,) => 0,
        (1, 3) => 2,
        (1, 4) => 3,
        (3, 4) => 6.0,
        (1, 2, 3) => 7,
        (1, 2, 4) => 8,
        (1, 3, 4) => 9,
        (1, 2, 3, 4) => 9,
        (3,) => 10_000,
    ]; threshold=8)

    @test flt isa Custom{Float64}
    @test dim(flt) == 3
    @test sort(flt[0], by=index) == [
        Simplex{0}(1, 0.0),
        Simplex{0}(2, 7.0),
        Simplex{0}(3, 2.0),
        Simplex{0}(4, 0.0),
    ]
    @test simplex(flt, Val(2), (1, 2, 4)) === Simplex{2}((4, 2, 1), 8.0)
    @test simplex(flt, Val(2), (4, 3, 1)) === nothing

    d0, d1, d2 = ripserer(flt, dim_max=2)
    @test d0 == [(0, 3), (0, Inf)]
    @test d1 == [(6, Inf)]
    @test d2 == []
end
