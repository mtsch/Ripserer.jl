using Random
using Ripserer
using Test

using Ripserer: chain_element_type, coefficient, index

using Ripserer: ReducedMatrix, record!, commit!, clear_buffer!
using Ripserer: WorkingChain, nonheap_push!, repair!

cofacet_type(::Type{<:A}) where {D,T,I,A<:Simplex{D,T,I}} = Simplex{D + 1,T,I}
facet_type(::Type{<:A}) where {D,T,I,A<:Simplex{D,T,I}} = Simplex{D - 1,T,I}

@testset "ReducedMatrix" begin
    for S in (Simplex{2,Int,Int}, Simplex{1,Int,Int32}), T in (Mod{3}, Rational{Int})
        C = cofacet_type(S)
        CE = chain_element_type(C, T)
        SE = chain_element_type(S, T)

        columns = C.([3, 10, 6, 8], 1)
        colelems = CE.(columns)

        fwd = Base.Order.Forward
        rev = Base.Order.Reverse

        @testset "ReducedMatrix with $S and $T" begin
            @testset "a fresh ReducedMatrix is empty" begin
                matrix = ReducedMatrix{C,SE}(fwd)
                @test length(matrix) == 0
                for col in Iterators.flatten((columns, colelems))
                    @test isempty(matrix[col])
                end
            end

            @testset "is empty after committing nothing" begin
                matrix = ReducedMatrix{C,SE}(fwd)
                commit!(matrix, columns[1], T(2))

                @test length(matrix) == 0
                for col in Iterators.flatten((columns, colelems))
                    @test isempty(matrix[col])
                end
            end

            @testset "is still empty after recording values" begin
                matrix = ReducedMatrix{C,SE}(rev)
                vals = [SE(S(1, 1)), SE(S(2, 1)), SE(S(3, 1)), SE(S(4, 1))]
                record!(matrix, vals, one(T))
                record!(matrix, S(5, 1))

                @test length(matrix) == 0
                for col in Iterators.flatten((columns, colelems))
                    @test isempty(matrix[col])
                end
            end

            @testset "commiting adds a column and changes length" begin
                matrix = ReducedMatrix{C,SE}(fwd)
                vals = [SE(S(1, 1)), SE(S(2, 1)), SE(S(3, 1)), SE(S(4, 1)), SE(S(-1, 1))]
                record!(matrix, vals, one(T))
                record!(matrix, S(4, 1))

                commit!(matrix, columns[1], T(2))
                clear_buffer!(matrix)
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

            @testset "clear_buffer! undos changes" begin
                matrix = ReducedMatrix{C,SE}(rev)
                vals = [
                    SE(S(1, 1)),
                    SE(S(2, 1)),
                    SE(S(3, 1)),
                    SE(S(4, 1)),
                    SE(S(4, 1)),
                    SE(S(-1, 1)),
                ]
                record!(matrix, vals, one(T))
                clear_buffer!(matrix)
                commit!(matrix, columns[2], T(2))

                record!(matrix, vals, one(T))
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
                matrix = ReducedMatrix{C,SE}(fwd)
                vals = [
                    SE(S(1, 1)),
                    SE(S(-2, 1)),
                    SE(S(4, 1)),
                    SE(S(-4, 1)),
                    SE(S(2, 1)),
                    SE(S(-1, 1)),
                ]
                record!(matrix, vals, one(T))
                commit!(matrix, columns[1], T(2))

                @test length(matrix) == 0
                for col in Iterators.flatten((columns, colelems))
                    @test isempty(matrix[col])
                end
            end

            @testset "committing multiple times creates multiple columns" begin
                matrix = ReducedMatrix{C,SE}(rev)
                vals_1 = [
                    SE(S(-1, 1)), SE(S(4, 1)), SE(S(-2, 1)), SE(S(-3, 1)), SE(S(-4, 1))
                ]
                record!(matrix, vals_1, -one(T))
                commit!(matrix, columns[1], T(2))
                clear_buffer!(matrix)

                vals_2 = [SE(S(4, 1)), SE(S(1, 1)), SE(S(5, 1)), SE(S(6, 1)), SE(S(-1, 1))]
                record!(matrix, vals_2, one(T))
                commit!(matrix, columns[2], T(-1))
                clear_buffer!(matrix)

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

@testset "WorkingChain" begin
    for S in (Simplex{2,Int,Int}, Simplex{1,Int,Int32}), T in (Mod{3}, Rational{Int})
        SE = chain_element_type(S, T)

        elements = SE.(S.([1, -7, 2, 3, 4, 7, 5, 6, -1], [1, 7, 1, 1, 4, 7, 5, 6, 1]))
        unq_elements = SE.(S.([7, 2, 3, 4, 5, 6], [7, 1, 1, 4, 5, 6]))

        fwd = Base.Order.Forward
        rev = Base.Order.Reverse

        @testset "a fresh WorkingChain is empty and pop! yields nothing" begin
            working_boundary = WorkingChain{SE}(fwd)

            @test isempty(working_boundary)
            @test pop!(working_boundary) ≡ nothing
            @test_throws BoundsError first(working_boundary)
        end

        @testset "pushing elements and popping finds the lowest simplex (cohomology)" begin
            working_boundary = WorkingChain{SE}(fwd)
            for e in elements
                push!(working_boundary, e)
            end

            @test pop!(working_boundary) == SE(S(3, 1))
        end

        @testset "the same happens for nonheap_push! with repair! (homology)" begin
            working_boundary = WorkingChain{SE}(rev)
            for e in elements
                nonheap_push!(working_boundary, e)
            end
            repair!(working_boundary)

            @test pop!(working_boundary) == SE(S(6, 6))
        end

        @testset "adding inverses removes elements" begin
            working_boundary = WorkingChain{SE}(fwd)
            for e in unq_elements
                nonheap_push!(working_boundary, e)
            end
            repair!(working_boundary)
            for e in unq_elements[1:2:end]
                push!(working_boundary, -e)
            end

            for e in sort(unq_elements[2:2:end]; order=fwd)
                @test pop!(working_boundary) == e
            end
            @test pop!(working_boundary) ≡ nothing
        end
    end
end
