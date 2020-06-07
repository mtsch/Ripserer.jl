using Ripserer

using Random

using Ripserer: chain_element_type, coefficient

using Ripserer: ReducedMatrix, record!, commit!, discard!
using Ripserer: WorkingBoundary, nonheap_push!, get_pivot!, repair!
using Ripserer: ReductionMatrix, simplex_type, simplex_element, face_element

@testset "ReducedMatrix" begin
    for S in (Simplex{2, Int, Int}, Cubelet{1, Int, Int}), T in (Mod{3}, Rational{Int})
        C = coface_type(S)
        CE = chain_element_type(C, T)
        SE = chain_element_type(S, T)

        columns = C.([3, 10, 6, 8], 1)
        colelems = CE.(columns)

        @testset "ReducedMatrix with simplex type $S and field type $T" begin
            @testset "a fresh ReducedMatrix is empty even if you commit nothing." begin
                matrix = ReducedMatrix{C, SE}(true)
                commit!(matrix, columns[1], T(2))

                @test length(matrix) == 0

                for col in Iterators.flatten((columns, colelems))
                    @test isempty(matrix[col])
                end
            end

            @testset "is still empty after recording values" begin
                matrix = ReducedMatrix{C, SE}(false)
                vals = [
                    SE(S(1, 1)),
                    SE(S(2, 1)),
                    SE(S(3, 1)),
                    SE(S(4, 1)),
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
                matrix = ReducedMatrix{C, SE}(true)
                vals = [
                    SE(S(1, 1)),
                    SE(S(2, 1)),
                    SE(S(3, 1)),
                    SE(S(4, 1)),
                    SE(S(4, 1)),
                    SE(S(-1, 1)),
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
                matrix = ReducedMatrix{C, SE}(false)
                vals = [
                    SE(S(1, 1)),
                    SE(S(2, 1)),
                    SE(S(3, 1)),
                    SE(S(4, 1)),
                    SE(S(4, 1)),
                    SE(S(-1, 1)),
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
                matrix = ReducedMatrix{C, SE}(true)
                vals = [
                    SE(S(1, 1)),
                    SE(S(-2, 1)),
                    SE(S(4, 1)),
                    SE(S(-4, 1)),
                    SE(S(2, 1)),
                    SE(S(-1, 1)),
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
                matrix = ReducedMatrix{C, SE}(false)
                vals_1 = [
                    SE(S(1, 1)),
                    SE(S(-4, 1)),
                    SE(S(2, 1)),
                    SE(S(3, 1)),
                    SE(S(4, 1)),
                ]
                for v in vals_1
                    record!(matrix, v)
                end
                commit!(matrix, columns[1], T(2))

                vals_1 = [
                    SE(S(4, 1)),
                    SE(S(1, 1)),
                    SE(S(5, 1)),
                    SE(S(6, 1)),
                    SE(S(-1, 1)),
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
end

@testset "WorkingBoundary" begin
    for S in (Simplex{2, Int, Int}, Cubelet{1, Int, Int}), T in (Mod{3}, Rational{Int})
        SE = chain_element_type(S, T)

        elements = SE.(S.([1, -7, 2, 3, 4, 7, 5, 6, -1], [1, 7, 1, 1, 4, 7, 5, 6, 1]))

        @testset "a fresh WorkingBoundary is empty and its pivot is nothing" begin
            working_boundary = WorkingBoundary{SE}(true)

            @test isempty(working_boundary)
            @test get_pivot!(working_boundary) ≡ nothing
            @test first(working_boundary) ≡ nothing
        end

        @testset "pushing elements and getting pivot finds the lowest simplex (cohomology)" begin
            working_boundary = WorkingBoundary{SE}(true)
            for e in elements
                push!(working_boundary, e)
            end

            @test get_pivot!(working_boundary) == SE(S(3, 1))
        end

        @testset "the same happens for nonheap_push! with repair! (homology)" begin
            working_boundary = WorkingBoundary{SE}(false)
            for e in elements
                nonheap_push!(working_boundary, e)
            end
            repair!(working_boundary)

            @test get_pivot!(working_boundary) == SE(S(6, 6))
        end

        @testset "getting pivot and adding its inverse to the boundary repeatedly" begin
            working_boundary = WorkingBoundary{SE}(true)
            for e in elements[1:3]
                nonheap_push!(working_boundary, e)
            end
            repair!(working_boundary)
            for e in elements[4:end]
                push!(working_boundary, e)
            end

            result = []
            while (p = get_pivot!(working_boundary)) ≢ nothing
                push!(working_boundary, -p)
                push!(result, p)
            end
            @test length(result) == 5
            @test allunique(result)
            @test issorted(result)
        end
    end
end

@testset "ReducedMatrix" begin
    for Co in (true, false),
        S in (Simplex{2, Int, Int}, Cubelet{1, Int, Int}),
        T in (Mod{3}, Rational{Int})

        columns_1 = (S.([1, 2, 3, 4, 5], [1, 2, 3, 4, 5]))
        columns_2 = (S.([6, 7, 8, 9, 10], [6, 7, 8, 9, 10]))
        flt = "anything really"
        co_str = Co ? "co" : ""

        SE = chain_element_type(S, T)
        F = Co ? coface_type(S) : face_type(S)
        FE = chain_element_type(F, T)

        @testset "$(co_str)homology, $S, $T" begin
            @testset "construction inferrence, parameters" begin
                @test begin
                    @inferred ReductionMatrix{Co, T}(flt, columns_1, columns_2)
                    true end

                matrix = ReductionMatrix{Co, T}(flt, columns_1, columns_2)

                @test simplex_type(matrix) == S
                @test simplex_element(matrix) == SE
                @test face_element(matrix) == FE
                @test dim(matrix) == (Co ? dim(S) : dim(S) - 1)
            end
        end
    end
end
