using Ripserer
using Test

using Ripserer: Chain, ChainElement, ReducedMatrix, bulk_insert!

for F in (Mod{2}, Mod{3}, Mod{257}, Rational{Int})
    @testset "ReducedMatrix with $F" begin
        I = Simplex{2,Int,Int}
        S = Simplex{3,Int,Int}

        columns = I.([3, 10, 6, 8], 1)
        elements = ChainElement{I,F}.(columns)

        @testset "a fresh ReducedMatrix is empty" begin
            @test begin
                @inferred ReducedMatrix{I,F,S}()
                @inferred ReducedMatrix{S,F,I}()
                true
            end

            matrix = ReducedMatrix{I,F,S}()
            @test length(matrix) == 0
            for col in Iterators.flatten((columns, elements))
                @test isempty(matrix[col])
            end
        end

        @testset "`setindex!` adds a column and changes length" begin
            matrix = ReducedMatrix{I,F,S}()
            vals1 = Chain{F,S}([S(1, 1), S(2, 1), S(3, 1), S(4, 1)])

            matrix[columns[1]] = vals1
            @test length(matrix) == 1
            @test matrix[elements[1]] == vals1

            for col in Iterators.flatten((columns[2:end], elements[2:end]))
                @test isempty(matrix[col])
            end

            vals2 = Chain{F,S}([S(1, 1), S(2, 1)])
            matrix[columns[2]] = vals2

            @test length(matrix) == 2
            @test matrix[elements[2]] == vals2

            for col in Iterators.flatten((columns[3:end], elements[3:end]))
                @test !haskey(matrix, col)
            end
        end

        @testset "`bulk_insert!`" begin
            matrix = ReducedMatrix{I,F,S}()
            bulk_insert!(matrix, ())
            @test length(matrix) == 0

            pairs = [
                (S(1, 1), I(2, 2)),
                (S(3, 3), I(4, 4)),
                (S(5, 5), I(6, 6)),
                (S(7, 7), I(8, 8))
            ]
            bulk_insert!(matrix, pairs)

            for (σ, τ) in pairs
                @test matrix[τ] == Chain{F,S}([σ])
            end
            @test length(matrix) == 4
        end
    end
end
