using RecipesBase: apply_recipe
using Ripserer: Barcode

# Hack to avoid having to import Plots.
RecipesBase.is_key_supported(::Symbol) = true

@testset "diagrams" begin
    @testset "PersistenceInterval" begin
        @testset "no cocycle" begin
            int1 = PersistenceInterval(1, 2)
            @test eltype(int1) == Union{Int, Infinity}
            b, d = int1
            @test b == birth(int1) == 1
            @test d == death(int1) == 2
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

            @test int1 < int2
            @test int1 < PersistenceInterval(2, 2)

            @test isnothing(cocycle(int1))

            @test sprint(print, int1) == "[1, 2)"
            @test sprint(print, int2) == "[1, ∞)"
            @test sprint((io, val) -> show(io, MIME"text/plain"(), val), int1) ==
                "PersistenceInterval{Int64}(1, 2)"
        end

        @testset "with cocycle" begin
            int1 = PersistenceInterval(2.0, 3.0, [1, 2, 3, 4])
            @test eltype(int1) == Union{Float64, Infinity}
            b, d = int1
            @test b == birth(int1) == 2.0
            @test d == death(int1) == 3.0
            @test int1 == (2.0, 3.0)
            @test_throws BoundsError int1[0]
            @test int1[1] == 2.0
            @test int1[2] == 3.0
            @test_throws BoundsError int1[3]

            int2 = PersistenceInterval(1.0, ∞, [1, 2])
            @test eltype(int2) == Union{Float64, Infinity}
            b, d = int2
            @test b == birth(int2) == 1.0
            @test d == death(int2) == ∞

            @test int1 > int2
            @test int1 > PersistenceInterval(2.0, 2.0, [1, 2])

            @test cocycle(int1) == [1, 2, 3, 4]
            @test cocycle(int2) == [1, 2]

            @test sprint(print, int1) == "[2.0, 3.0)"
            @test sprint(print, int2) == "[1.0, ∞)"
            @test sprint((io, val) -> show(io, MIME"text/plain"(), val), int1) ==
                """
                PersistenceInterval{Float64}(2.0, 3.0) with cocycle:
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
            int2 = PersistenceInterval(1, 2)
            int3 = PersistenceInterval(3, 4)

            diag = PersistenceDiagram(1, [int1, int2, int3])
            @test dim(diag) == 1
            @test diag ==  [int2, int3, int1]
            @test diag[1] == int2
            @test diag[2] == int3
            @test diag[3] == int1
            @test_throws BoundsError diag[4]
            @test sprint(print, diag) ==
                "3-element 1-dimensional PersistenceDiagram"
            @test sprint((io, val) -> show(io, MIME"text/plain"(), val), diag) ==
                """
                3-element 1-dimensional PersistenceDiagram:
                 [1, 2)
                 [3, 4)
                 [3, ∞)"""
        end
        @testset "plots recipes" begin
            int1 = PersistenceInterval(3, ∞)
            int2 = PersistenceInterval(1, 2)
            int3 = PersistenceInterval(3, 4)

            diag1 = PersistenceDiagram(1, [int1, int2, int3])
            diag2 = PersistenceDiagram(1, [int2, int3])

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
end
