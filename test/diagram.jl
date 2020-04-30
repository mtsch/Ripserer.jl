using Ripserer

using Compat
using RecipesBase
using RecipesBase: apply_recipe
using Ripserer: Barcode

# Hack to avoid having to import Plots.
RecipesBase.is_key_supported(::Symbol) = true

@testset "PersistenceInterval" begin
    @testset "no representative" begin
        int1 = PersistenceInterval(1, 2)
        @test eltype(int1) == Union{Int, Infinity}
        b, d = int1
        @test b == birth(int1) == 1
        @test d == death(int1) == 2
        @test persistence(int1) == 1
        @test int1 == (1, 2)
        @test_throws BoundsError int1[0]
        @test int1[1] == 1
        @test int1[2] == 2
        @test_throws BoundsError int1[3]

        int2 = PersistenceInterval(1, ∞)
        @test eltype(int2) == Union{Int, Infinity}
        b, d = int2
        @test b == birth(int2) == 1
        @test d == death(int2) == ∞
        @test persistence(int2) == ∞

        @test int1 < int2
        @test int1 < PersistenceInterval(2, 2)

        @test isnothing(representative(int1))

        @test sprint(print, int1) == "[1, 2)"
        @test sprint(print, int2) == "[1, ∞)"
        @test sprint((io, val) -> show(io, MIME"text/plain"(), val), int1) ==
            "PersistenceInterval{Int64}(1, 2)"
    end
    @testset "with representative" begin
        int1 = PersistenceInterval(2.0, 3.0, [1, 2, 3, 4])
        @test eltype(int1) == Union{Float64, Infinity}
        b, d = int1
        @test b == birth(int1) == 2.0
        @test d == death(int1) == 3.0
        @test persistence(int1) == 1.0
        @test int1 == (2.0, 3.0)
        @test_throws BoundsError int1[0]
        @test int1[1] == 2.0
        @test int1[2] == 3.0
        @test_throws BoundsError int1[3]

        int2 = PersistenceInterval(1.0, ∞, [1, 2])
        @test eltype(int2) == Union{Float64, Infinity}
        b, d = int2
        @test b == birth(int2) == 1.0
        @test d == death(int2) == Inf
        @test persistence(int2) == Inf

        @test int1 > int2
        @test int1 > PersistenceInterval(2.0, 2.0, [1, 2])

        @test representative(int1) == [1, 2, 3, 4]
        @test representative(int2) == [1, 2]

        @test sprint(print, int1) == "[2.0, 3.0)"
        @test sprint(print, int2) == "[1.0, ∞)"
        @test sprint((io, val) -> show(io, MIME"text/plain"(), val), int1) ==
            """
                PersistenceInterval{Float64}(2.0, 3.0) with representative:
                4-element Array{Int64,1}:
                 1
                 2
                 3
                 4"""
    end
end

@testset "PersistenceDiagram" begin
    @testset "basics" begin
        int1 = PersistenceInterval(3, ∞)
        int2 = PersistenceInterval(1, 3)
        int3 = PersistenceInterval(3, 4)

        diag = PersistenceDiagram(1, [int1, int2, int3])
        @test dim(diag) == 1
        sort!(diag)
        @test diag == [int2, int3, int1]
        @test diag == PersistenceDiagram(1, [(1, 3), (3, 4), (3, ∞)])
        @test diag[1] == int2
        @test diag[2] == int3
        @test diag[3] == int1
        @test_throws BoundsError diag[4]
        @test sprint(print, diag) ==
            "3-element 1-dimensional PersistenceDiagram"
        @test sprint((io, val) -> show(io, MIME"text/plain"(), val), diag) ==
            """
                3-element 1-dimensional PersistenceDiagram:
                 [1, 3)
                 [3, 4)
                 [3, ∞)"""

        @test copy(diag) == diag
        @test similar(diag) isa typeof(diag)
        @test similar(diag, (Base.OneTo(2),)) isa typeof(diag)
        @test sort(diag) isa typeof(diag)
        @test filter(isfinite, diag) isa typeof(diag)
    end
    @testset "Plots recipes" begin
        # Idea: apply recipe and check the number of series on plots.
        # Not a perfect way to test, but at least it makes sure the infinity is drawn only
        # when it needs to be and that there are no errors.
        int1 = PersistenceInterval(3, ∞)
        int2 = PersistenceInterval(1, 2)
        int3 = PersistenceInterval(3, 4)

        diag1 = PersistenceDiagram(1, [int1, int2, int3])
        diag2 = PersistenceDiagram(2, [(1, 2), (3, 4)])

        d = Dict{Symbol, Any}()
        @test length(apply_recipe(d, diag1)) == 3
        d = Dict{Symbol, Any}()
        @test length(apply_recipe(d, diag2)) == 2
        d = Dict{Symbol, Any}()
        @test length(apply_recipe(d, [diag1, diag2])) == 4

        d = Dict{Symbol, Any}()
        @test length(apply_recipe(d, Barcode((diag1,)))) == 2
        d = Dict{Symbol, Any}()
        @test length(apply_recipe(d, Barcode((diag2,)))) == 1
        d = Dict{Symbol, Any}()
        @test length(apply_recipe(d, Barcode(([diag1, diag2],)))) == 3

        d = Dict{Symbol, Any}(:infinity => 10)
        @test length(apply_recipe(d, diag1)) == 3
        d = Dict{Symbol, Any}(:infinity => 10)
        @test length(apply_recipe(d, diag2)) == 2
        d = Dict{Symbol, Any}(:infinity => 10)
        @test length(apply_recipe(d, [diag1, diag2])) == 4

        d = Dict{Symbol, Any}(:infinity => 10)
        @test length(apply_recipe(d, Barcode((diag1,)))) == 2
        d = Dict{Symbol, Any}(:infinity => 10)
        @test length(apply_recipe(d, Barcode((diag2,)))) == 1
        d = Dict{Symbol, Any}(:infinity => 10)
        @test length(apply_recipe(d, Barcode(([diag1, diag2],)))) == 3

        d = Dict{Symbol, Any}(:infinity => 10)
        @test_throws ArgumentError apply_recipe(d, Barcode((1,)))
        d = Dict{Symbol, Any}()
        @test_throws ArgumentError apply_recipe(d, Barcode((diag1,2)))
    end
end
