using PersistenceDiagrams
using Ripserer
using RecipesBase
using Test

using Ripserer: plottable, ChainElement

using RecipesBase: apply_recipe

# Hack to avoid having to import Plots.
RecipesBase.is_key_supported(::Symbol) = true

"""
Idea: apply recipe and check the number of series on plots.
Not a perfect way to test, but at least it makes sure all points are plotted and checks that
there are no errors.
"""
series(args...; kwargs...) =
    apply_recipe(Dict{Symbol, Any}(kwargs), args...)

@testset "simplex plots" begin
    @testset "plottable" begin
        @test isequal(
            plottable(Simplex{0}(3, 1), 1:10),
            ([3],
             [:seriestype => :scatter], 0),
        )
        @test isequal(
            plottable(Simplex{1}(1, 1),  1:10),
            ([2, 1, NaN],
             [:seriestype => :path], 1),
        )
        @test isequal(
            plottable(Simplex{2}(3, 1), 1:10),
            ([4, 3, 1, 4, NaN],
             [:seriestype => :path], 2),
        )
        @test isequal(
            plottable(Simplex{3}(1, 1), 1:10),
            ([4, 3, 2, 4, 4, 3, 1, 4, 4, 2, 1, 4, 3, 2, 1, 3, NaN],
             [:seriestype => :path], 3),
        )
    end
    @testset "plot" begin
        # No idea how to test this better. At least it checks there are no erros.
        data = collect(1:100)
        for dim in 0:3
            sx = Simplex{dim}(1, 1)
            fsx = Simplex{dim + 1}(1, 1)
            @test length(series(sx, data)) == 1
            @test length(series([sx], data)) == 1
            @test length(series([sx], data)) == 1
            @test length(series(RepresentativeInterval(
                1.0, 1.0, sx, fsx, [ChainElement{typeof(sx), typeof(1//1)}(sx, 1//1)]
            ), data)) == 1
            @test_throws ErrorException series(PersistenceInterval(1.0, 1.0), data)
        end
    end
end
