using Ripserer

using Ripserer: ReducedMatrix, record!, commit!, discard!
using Ripserer: ChainElement, coefficient

for T in (Mod{3}, Rational{Int})
    CE = ChainElement{Simplex{3, Int, Int}, T}
    SE = ChainElement{Simplex{2, Int, Int}, T}

    columns = Simplex{3}.([3, 10, 6, 8], 1)
    colelems = CE.(columns)

    @testset "Reduction matrix with $T" begin
        @testset "a fresh ReducedMatrix is empty even if you commit nothing." begin
            matrix = ReducedMatrix{Simplex{3, Int, Int}, SE}(false)
            commit!(matrix, columns[1], T(2))

            @test length(matrix) == 0

            for col in Iterators.flatten((columns, colelems))
                @test isempty(matrix[col])
            end
        end

        @testset "a ReducedMatrix is still empty after recording values" begin
            matrix = ReducedMatrix{Simplex{3, Int, Int}, SE}(true)
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
            matrix = ReducedMatrix{Simplex{3, Int, Int}, SE}(false)
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

        @testset "discard! undos changes" begin
            matrix = ReducedMatrix{Simplex{3, Int, Int}, SE}(true)
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
            discard!(matrix)
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
            matrix = ReducedMatrix{Simplex{3, Int, Int}, SE}(false)
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

        @testset "committing multiple times creates multiple columns" begin
            matrix = ReducedMatrix{Simplex{3, Int, Int}, SE}(true)
            vals_1 = [
                SE(Simplex{2}(1, 1)),
                SE(Simplex{2}(-4, 1)),
                SE(Simplex{2}(2, 1)),
                SE(Simplex{2}(3, 1)),
                SE(Simplex{2}(4, 1)),
            ]
            for v in vals_1
                record!(matrix, v)
            end
            commit!(matrix, columns[1], T(2))

            vals_1 = [
                SE(Simplex{2}(4, 1)),
                SE(Simplex{2}(1, 1)),
                SE(Simplex{2}(5, 1)),
                SE(Simplex{2}(6, 1)),
                SE(Simplex{2}(-1, 1)),
            ]
            for v in vals_1
                record!(matrix, v)
            end
            commit!(matrix, columns[2], T(-1))

            @test length(matrix) == 2

            @test index.(simplex.(collect(matrix[columns[1]]))) == [1, 2, 3]
            @test index.(simplex.(collect(matrix[columns[2]]))) == [4, 5, 6]

            for col in Iterators.flatten((columns[3:end], colelems[3:end]))
                @test isempty(matrix[col])
            end
        end
    end
end
