using Ripserer: DiameterSimplexHeap, CurrentColumn, initialize!, pop_pivot!

@testset "columns" begin
    @testset "DiameterSimplexHeap" begin
        column = DiameterSimplexHeap{3}()
        @test isempty(column)

        push!(column, DiameterSimplex{3}(3.0, 1, 1))
        push!(column, DiameterSimplex{3}(3.0, 2, 1))
        push!(column, DiameterSimplex{3}(1.0, 3, 1))
        push!(column, DiameterSimplex{3}(2.0, 4, 1))
        push!(column, DiameterSimplex{3}(2.0, 4, 1))
        @test length(column) == 5
        @test !isempty(column)

        @test top(column) == DiameterSimplex{3}(1.0, 3, 1)
        @test pop!(column) == DiameterSimplex{3}(1.0, 3, 1)
        @test length(column) == 4
        @test pop!(column) == DiameterSimplex{3}(2.0, 4, 1)
        @test pop!(column) == DiameterSimplex{3}(2.0, 4, 1)
        @test pop!(column) == DiameterSimplex{3}(3.0, 2, 1)
        @test pop!(column) == DiameterSimplex{3}(3.0, 1, 1)
        @test isempty(column)

        @test_throws MethodError push!(column, DiameterSimplex{2}(3.0, 1, 1))
    end

    @testset "initialize!" begin
        @testset "star" begin
            col = CurrentColumn{5}()
            dist = sparse([0 1 0 0;
                           1 0 2 3;
                           0 2 0 0;
                           0 3 0 0])

            initialize!(col, DiameterSimplex{5}(1, [2], 1), 0, dist, binomial)
            @test length(col.heap) == 3
            @test pop!(col.heap) == DiameterSimplex{5}(1, 1, 1)
            @test pop!(col.heap) == DiameterSimplex{5}(2, 3, 1)
            @test pop!(col.heap) == DiameterSimplex{5}(3, 5, 1)
        end

        @testset "line cofaces" begin
            col = CurrentColumn{2}()
            dist = sparse(Float64[0 1 3 4 5 0;
                                  1 0 3 4 5 1;
                                  3 3 0 0 0 1;
                                  4 4 0 0 0 1;
                                  5 5 0 0 0 1;
                                  0 1 1 1 1 1])
            initialize!(col, DiameterSimplex{2}(1.0, [2, 1], 1), 1, dist, binomial)
            @test length(col.heap) == 3
            @test pop!(col.heap) == DiameterSimplex{2}(3.0, [3, 2, 1], 1)
            @test pop!(col.heap) == DiameterSimplex{2}(4.0, [4, 2, 1], 1)
            @test pop!(col.heap) == DiameterSimplex{2}(5.0, [5, 2, 1], 1)
        end

        @testset "full graph" begin
            dist = ones(100, 100)
            for i in 1:size(dist, 1)
                dist[i, i] = 0
            end
            dist_sp = sparse(dist)

            for dim in 1:5
                col = CurrentColumn{3}()
                initialize!(col, DiameterSimplex{3}(1.0, 1, 1), dim, dist, binomial)
                @test length(col.heap) == 100 - dim - 1

                col = CurrentColumn{3}()
                initialize!(col, DiameterSimplex{3}(1.0, 1, 1), dim, dist_sp, binomial)
                @test length(col.heap) == 100 - dim - 1
            end
        end
    end

    @testset "pop_pivot!" begin
        @testset "single element" begin
            col = CurrentColumn{2}()
            push!(col.heap, DiameterSimplex{2}(2.0, 1, 1))
            push!(col.heap, DiameterSimplex{2}(2.0, 1, 1))
            push!(col.heap, DiameterSimplex{2}(2.0, 1, 1))
            push!(col.heap, DiameterSimplex{2}(2.0, 1, 1))
            push!(col.heap, DiameterSimplex{2}(2.0, 1, 1))

            @test pop_pivot!(col) == DiameterSimplex{2}(2.0, 1, 1)
            @test isempty(col.heap)

            col = CurrentColumn{3}()
            push!(col.heap, DiameterSimplex{3}(2.0, 1, 1))
            push!(col.heap, DiameterSimplex{3}(2.0, 1, 1))
            push!(col.heap, DiameterSimplex{3}(2.0, 1, 1))

            @test isnothing(pop_pivot!(col))
            @test isempty(col.heap)
        end

        @testset "multiple" begin
            col = CurrentColumn{5}()
            push!(col.heap, DiameterSimplex{5}(1.0, 2, 3))
            push!(col.heap, DiameterSimplex{5}(2.0, 3, 4))
            push!(col.heap, DiameterSimplex{5}(1.0, 2, 2))
            push!(col.heap, DiameterSimplex{5}(3.0, 1, 2))
            push!(col.heap, DiameterSimplex{5}(2.0, 3, 1))
            push!(col.heap, DiameterSimplex{5}(4.0, 4, 4))
            push!(col.heap, DiameterSimplex{5}(4.0, 4, 4))
            push!(col.heap, DiameterSimplex{5}(4.0, 4, 4))

            @test pop_pivot!(col) == DiameterSimplex{5}(3.0, 1, 2)
            @test pop_pivot!(col) == DiameterSimplex{5}(4.0, 4, 2)
            @test isnothing(pop_pivot!(col))
        end
    end
end
