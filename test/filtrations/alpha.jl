using Test
using Ripserer
using Suppressor
using StaticArrays

using Ripserer: circumcenter_radius2

include("interfacetest.jl")

@static if Sys.iswindows()
    @test_broken begin
        Alpha(points)
        true
    end
else
    @testset "circumcenter_radius2" begin
        c1, r2 = circumcenter_radius2([SA[1.0, 2.0, 3.0, 4.0]])
        @test c1 == [1.0, 2.0, 3.0, 4.0]
        @test r2 == 0

        c2, r2 = circumcenter_radius2([SA[1, 1, 1], SA[-1, -1, -1]])
        @test c2 == [0, 0, 0]
        @test r2 == 3.0

        c3, r3 = circumcenter_radius2([SA[1, 0], SA[-1, 0], SA[0, √3]])
        @test c3 ≈ [0, 1 / √3]
        @test r3 ≈ (2 / √3)^2

        c4, r4 = circumcenter_radius2([SA[1, 0], SA[0, 0], SA[-1, 0]])
        c4 = [NaN, NaN]
        r4 = Inf
    end

    @testset "Compare to Cechmate" begin
        # This test compares the results to cechmate: https://github.com/scikit-tda/cechmate/
        points = [(t * sinpi(t), t * cospi(t)) for t in range(0, 2; length=20)]

        d0_cechmate_deaths =
            2 .* [
                0.05263158,
                0.05805554,
                0.06761044,
                0.07982644,
                0.09366810,
                0.10851510,
                0.12400688,
                0.13992943,
                0.15615105,
                0.17258742,
                0.18918257,
                0.20589813,
                0.22270697,
                0.23958947,
                0.25653108,
                0.27352082,
                0.29055025,
                0.30761277,
                0.32470317,
                Inf,
            ]
        d1_cechmate_births = 2 .* [0.83075319, 0.86707034, 0.91920916, 0.89414800]
        d1_cechmate_deaths = 2 .* [0.83143060, 0.86803555, 0.92415392, 0.89985412]

        d0, d1, d2 = ripserer(Alpha(points); dim_max=2)
        @test all(iszero, birth.(d0))
        @test death.(d0) ≈ d0_cechmate_deaths atol = 1e-7
        @test birth.(d1) ≈ d1_cechmate_births atol = 1e-7
        @test death.(d1) ≈ d1_cechmate_deaths atol = 1e-7
        @test isempty(d2)
    end

    @testset "Interface" begin
        points = [(t * sinpi(t), t * cospi(t)) for t in range(0, 2; length=20)]
        test_filtration(Alpha, points)
    end
end
