using Ripserer

using Ripserer: ReductionMatrix, record!, commit!, undo!
using Ripserer: ChainElement, coefficient

for T in (Mod{3}, Rational{Int})
    CE = ChainElement{Simplex{3, Int, Int}, T}
    SE = ChainElement{Simplex{2, Int, Int}, T}

    columns = Simplex{3}.([3, 10, 6, 8], 1)
    colelems = CE.(columns)

    @testset "Reduction matrix with $T" begin
        @testset "a fresh ReductionMatrix is empty even if you commit nothing." begin
            matrix = ReductionMatrix{Simplex{3, Int, Int}, SE}(true)
            commit!(matrix, columns[1], T(2))

            @test length(matrix) == 0

            for col in Iterators.flatten((columns, colelems))
                @test isempty(matrix[col])
            end
        end

        @testset "a ReductionMatrix is still empty after recording values" begin
            matrix = ReductionMatrix{Simplex{3, Int, Int}, SE}(false)
            vals = [
                SE(Simplex{2}(1, 1)),
                SE(Simplex{2}(2, 1)),
                SE(Simplex{2}(3, 1)),
                SE(Simplex{2}(4, 1)),
            ]
            for v in vals
                record!(matrix, v)
            end

            @test length(matrix) == 0

            for col in Iterators.flatten((columns, colelems))
                @test isempty(matrix[col])
            end
        end

        @testset "commiting adds a column and changes length" begin
            matrix = ReductionMatrix{Simplex{3, Int, Int}, SE}(false)
            vals = [
                SE(Simplex{2}(1, 1)),
                SE(Simplex{2}(2, 1)),
                SE(Simplex{2}(3, 1)),
                SE(Simplex{2}(4, 1)),
                SE(Simplex{2}(4, 1)),
                SE(Simplex{2}(-1, 1)),
            ]
            for v in vals
                record!(matrix, v)
            end

            commit!(matrix, columns[1], T(2))
            commit!(matrix, columns[2], T(2))

            @test length(matrix) == 1
            @test index.(simplex.(collect(matrix[columns[1]]))) == [4, 3, 2]
            @test index.(simplex.(collect(matrix[colelems[1]]))) == [4, 3, 2]

            @test coefficient.(collect(matrix[columns[1]])) == [T(4), T(2), T(2)]
            @test coefficient.(collect(matrix[colelems[1]])) == [T(4), T(2), T(2)]

            for col in Iterators.flatten((columns[2:end], colelems[2:end]))
                @test isempty(matrix[col])
            end
        end

        @testset "undo! undos changes" begin
            matrix = ReductionMatrix{Simplex{3, Int, Int}, SE}(true)
            vals = [
                SE(Simplex{2}(1, 1)),
                SE(Simplex{2}(2, 1)),
                SE(Simplex{2}(3, 1)),
                SE(Simplex{2}(4, 1)),
                SE(Simplex{2}(4, 1)),
                SE(Simplex{2}(-1, 1)),
            ]
            for v in vals
                record!(matrix, v)
            end
            undo!(matrix)
            commit!(matrix, columns[2], T(2))

            for v in vals
                record!(matrix, v)
            end
            commit!(matrix, columns[1], T(2))

            @test length(matrix) == 1
            @test index.(simplex.(collect(matrix[columns[1]]))) == [2, 3, 4]
            @test index.(simplex.(collect(matrix[colelems[1]]))) == [2, 3, 4]

            @test coefficient.(collect(matrix[columns[1]])) == [T(2), T(2), T(4)]
            @test coefficient.(collect(matrix[colelems[1]])) == [T(2), T(2), T(4)]

            for col in Iterators.flatten((columns[2:end], colelems[2:end]))
                @test isempty(matrix[col])
            end
        end

        @testset "committing values that sum to 0 does not create a column" begin
            matrix = ReductionMatrix{Simplex{3, Int, Int}, SE}(true)
            vals = [
                SE(Simplex{2}(1, 1)),
                SE(Simplex{2}(-2, 1)),
                SE(Simplex{2}(4, 1)),
                SE(Simplex{2}(-4, 1)),
                SE(Simplex{2}(2, 1)),
                SE(Simplex{2}(-1, 1)),
            ]
            for v in vals
                record!(matrix, v)
            end
            commit!(matrix, columns[1], T(2))

            @test length(matrix) == 0
            for col in Iterators.flatten((columns, colelems))
                @test isempty(matrix[col])
            end
        end
    end
end
