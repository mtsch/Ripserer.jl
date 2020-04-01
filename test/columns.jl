using Ripserer: Column, coboundary, pop_pivot!,
    compute_0_dim_pairs!,
    ReductionMatrices, compute_pairs!

@testset "columns" begin
    @testset "a Column is a heap." begin
        column = Column{3, Float64}()
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
        @test_throws MethodError push!(column, DiameterSimplex{3}(3, 1, 1))
    end

    @testset "coboundary" begin
        @testset "start" begin
            st = ReductionState(sparse([0 1 0 0;
                                        1 0 2 3;
                                        0 2 0 0;
                                        0 3 0 0]), 0, 5)
            cb_set = Set{DiameterSimplex{5, Int}}()
            for sx in coboundary(st, DiameterSimplex{5}(1, 2, 1))
                push!(cb_set, sx)
            end
            @test cb_set == Set([DiameterSimplex{5}(1, 1, 1),
                                 DiameterSimplex{5}(2, 3, 1),
                                 DiameterSimplex{5}(3, 5, 1)])
        end

        @testset "line cofaces" begin
            st = ReductionState(sparse(Float64[0 1 3 4 5 0;
                                               1 0 3 4 5 1;
                                               3 3 0 0 0 1;
                                               4 4 0 0 0 1;
                                               5 5 0 0 0 1;
                                               0 1 1 1 1 0]), 1, 2)
            cb_set = Set{DiameterSimplex{2, Float64}}()
            for sx in coboundary(st, DiameterSimplex{2}(st, 1.0, (2, 1), 1))
                push!(cb_set, sx)
            end
            @test length(cb_set) == 3
            @test cb_set == Set([DiameterSimplex{2}(st, 3.0, [3, 2, 1], 1),
                                 DiameterSimplex{2}(st, 4.0, [4, 2, 1], 1),
                                 DiameterSimplex{2}(st, 5.0, [5, 2, 1], 1)])
        end

        @testset "full graph" begin
            dist = ones(100, 100)
            for i in 1:size(dist, 1)
                dist[i, i] = 0
            end
            st = ReductionState(sparse(dist), 5, 3)

            for d in 1:5
                st.dim[] = d
                cob = DiameterSimplex{3, Float64}[]
                for sx in coboundary(st, DiameterSimplex{3}(1.0, 10, 1))
                    push!(cob, sx)
                end
                @test length(cob) == 100 - d - 1
            end
        end
    end

    #=
    @testset "initialize!" begin
        @testset "star" begin
            col = Column{5, Int}()
            dist = sparse([0 1 0 0;
                           1 0 2 3;
                           0 2 0 0;
                           0 3 0 0])

            initialize!(col, DiameterSimplex{5}(1, [2], 1), 0, dist, binomial)
            @test length(col) == 3
            @test pop!(col) == DiameterSimplex{5}(1, 1, 1)
            @test pop!(col) == DiameterSimplex{5}(2, 3, 1)
            @test pop!(col) == DiameterSimplex{5}(3, 5, 1)
        end

        @testset "line cofaces" begin
            col = Column{2, Float64}()
            dist = sparse(Float64[0 1 3 4 5 0;
                                  1 0 3 4 5 1;
                                  3 3 0 0 0 1;
                                  4 4 0 0 0 1;
                                  5 5 0 0 0 1;
                                  0 1 1 1 1 1])
            initialize!(col, DiameterSimplex{2}(1.0, [2, 1], 1), 1, dist, binomial)
            @test length(col) == 3
            @test pop!(col) == DiameterSimplex{2}(3.0, [3, 2, 1], 1)
            @test pop!(col) == DiameterSimplex{2}(4.0, [4, 2, 1], 1)
            @test pop!(col) == DiameterSimplex{2}(5.0, [5, 2, 1], 1)
        end

        @testset "full graph" begin
            dist = ones(100, 100)
            for i in 1:size(dist, 1)
                dist[i, i] = 0
            end
            dist_sp = sparse(dist)

            for dim in 1:5
                col = Column{3, Float64}()
                initialize!(col, DiameterSimplex{3}(1.0, 1, 1), dim, dist, binomial)
                @test length(col) == 100 - dim - 1

                col = Column{3, Float64}()
                initialize!(col, DiameterSimplex{3}(1.0, 1, 1), dim, dist_sp, binomial)
                @test length(col) == 100 - dim - 1
            end
        end
    end
    =#

    @testset "pop_pivot!" begin
        @testset "single element" begin
            col = Column{2, Float64}()
            push!(col, DiameterSimplex{2}(2.0, 1, 1))
            push!(col, DiameterSimplex{2}(2.0, 1, 1))
            push!(col, DiameterSimplex{2}(2.0, 1, 1))
            push!(col, DiameterSimplex{2}(2.0, 1, 1))
            push!(col, DiameterSimplex{2}(2.0, 1, 1))

            @test pop_pivot!(col) == DiameterSimplex{2}(2.0, 1, 1)
            @test isempty(col)

            col = Column{3, Float64}()
            push!(col, DiameterSimplex{3}(2.0, 1, 1))
            push!(col, DiameterSimplex{3}(2.0, 1, 1))
            push!(col, DiameterSimplex{3}(2.0, 1, 1))

            @test isnothing(pop_pivot!(col))
            @test isempty(col)
        end

        @testset "multiple" begin
            col = Column{5, Float64}()
            push!(col, DiameterSimplex{5}(1.0, 2, 3))
            push!(col, DiameterSimplex{5}(2.0, 3, 4))
            push!(col, DiameterSimplex{5}(1.0, 2, 2))
            push!(col, DiameterSimplex{5}(3.0, 1, 2))
            push!(col, DiameterSimplex{5}(2.0, 3, 1))
            push!(col, DiameterSimplex{5}(4.0, 4, 4))
            push!(col, DiameterSimplex{5}(4.0, 4, 4))
            push!(col, DiameterSimplex{5}(4.0, 4, 4))

            @test pop_pivot!(col) == DiameterSimplex{5}(3.0, 1, 2)
            @test pop_pivot!(col) == DiameterSimplex{5}(4.0, 4, 2)
            @test isnothing(pop_pivot!(col))
        end
    end

    @testset "compute_0_dim_pairs!" begin
        @testset "dense Int" begin
            dist = [0 1 2;
                    1 0 3;
                    2 3 0]
            st = ReductionState(dist, 1, 2)
            critical_edges = DiameterSimplex{2, Int64}[]
            res = compute_0_dim_pairs!(st, critical_edges)

            @test res == [(0, 1),
                          (0, 2),
                          (0, typemax(Int))]
            @test critical_edges == [DiameterSimplex{2}(3, 3, 1)]
        end

        @testset "sparse Float64" begin
            dist = sparse(Float64[0 2 0 0 5 0;
                                  2 0 4 6 0 0;
                                  0 4 0 3 0 0;
                                  0 6 3 0 1 0;
                                  5 0 0 1 0 0;
                                  0 0 0 0 0 0])
            st = ReductionState(dist, 1, 3)
            critical_edges = DiameterSimplex{3, Float64}[]
            res = compute_0_dim_pairs!(st, critical_edges)

            @test res == [(0.0, 1.0),
                          (0.0, 2.0),
                          (0.0, 3.0),
                          (0.0, 4.0),
                          (0.0, Inf),
                          (0.0, Inf)]
            @test critical_edges == [DiameterSimplex{3}(6.0, 5, 1),
                                     DiameterSimplex{3}(5.0, 7, 1)]
        end
    end

    @testset "compute_pairs!" begin
        # 1 -- 2
        # |    |
        # 4 -- 3
        dist = Float64[0 1 2 1;
                       1 0 1 2;
                       2 1 0 1;
                       1 2 1 0]
        st = ReductionState(dist, 1, 2)
        columns = DiameterSimplex{2, Float64}[]
        compute_0_dim_pairs!(st, columns)

        rm = ReductionMatrices(st)
        res = compute_pairs!(rm, columns, 1)

        @test res == [(1.0, 2.0)]
    end

    @test "ripser" begin
        dist = rand_dist_matrix(100)
        st = ReductionState(dist, 1, 2)
        columns = DiameterSimplex{2, Float64}[]
        compute_0_dim_pairs!(st, columns)

        rmx = ReductionMatrices(st)
        res = compute_pairs!(rmx, columns, 1)
    end
end
