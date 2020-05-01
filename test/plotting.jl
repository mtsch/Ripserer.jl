using Ripserer

using RecipesBase
using RecipesBase: apply_recipe
using Ripserer: Barcode, plottable

# Hack to avoid having to import Plots.
RecipesBase.is_key_supported(::Symbol) = true

"""
Idea: apply recipe and check the number of series on plots.
Not a perfect way to test, but at least it makes sure the infinity is drawn only when it
needs to be and that there are no errors.
"""
series(args...; kwargs...) =
    apply_recipe(Dict{Symbol, Any}(kwargs), args...)

@testset "diagram plot, barcode" begin
    int1 = PersistenceInterval(3, ∞)
    int2 = PersistenceInterval(1, 2)
    int3 = PersistenceInterval(3, 4)

    diag1 = PersistenceDiagram(1, [int1, int2, int3])
    diag2 = PersistenceDiagram(2, [(1, 2), (3, 4)])

    # inf, x=y, points, ()
    @test length(series(diag1)) == 4
    # x=y, points, ()
    @test length(series(diag2)) == 3
    # inf, x=y, points 2×, ()
    @test length(series([diag1, diag2])) == 5

    # inf, lines, ()
    @test length(series(Barcode((diag1,)))) == 3
    # lines, ()
    @test length(series(Barcode((diag2,)))) == 2
    # inf, lines 2×, ()
    @test length(series(Barcode(([diag1, diag2],)))) == 4

    @test length(series(diag1, infinity=10)) == 4
    @test length(series(diag2, infinity=10)) == 3
    @test length(series([diag1, diag2], infinity=10)) == 5

    @test length(series(Barcode((diag1,)), infinity=10)) == 3
    @test length(series(Barcode((diag2,)), infinity=10)) == 2
    @test length(series(Barcode(([diag1, diag2],)), infinity=10)) == 4

    @test_throws MethodError series(Barcode((1,)))
    @test_throws MethodError series(Barcode((diag1,2)))
end

@testset "simplex plots" begin
    @testset "plottable" begin
        @test isequal(
            plottable(Simplex{0, 2}(1, 2, 3), 1:10),
            ([2], [:seriestype => :scatter]),
        )
        @test isequal(
            plottable(Simplex{1, 2}(1, 1, 3),  1:10),
            ([2, 1, NaN], [:seriestype => :path]),
        )
        @test isequal(
            plottable(Simplex{2, 2}(1, 3, 3), 1:10),
            ([4, 3, 1, 4, NaN], [:seriestype => :path]),
        )
        @test isequal(
            plottable(Simplex{3, 2}(1, 1, 3), 1:10),
            ([4, 3, 2, 4, 4, 3, 1, 4, 4, 2, 1, 4, 3, 2, 1, 3, NaN], [:seriestype => :path]),
        )
    end
    @testset "plot" begin
        # no idea how to test this better
        data = collect(1:100)
        for dim in 0:3
            sx = Simplex{dim, 2}(1, 1, 1)
            @test length(series(sx, data)) == 1
            @test length(series([sx], data)) == 1
            @test length(series([sx], data)) == 1
        end
    end
end
