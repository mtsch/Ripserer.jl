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
            st = ReductionState{5}(sparse([0 1 0 0;
                                           1 0 2 3;
                                           0 2 0 0;
                                           0 3 0 0]), 0)
            cb_set = Set{DiameterSimplex{5, Int}}()
            for sx in coboundary(st, DiameterSimplex{5}(1, 2, 1), 0)
                push!(cb_set, sx)
            end
            @test cb_set == Set([DiameterSimplex{5}(1, 1, 4),
                                 DiameterSimplex{5}(2, 3, 1),
                                 DiameterSimplex{5}(3, 5, 1)])
        end

        @testset "line cofaces" begin
            st = ReductionState{2}(sparse(Float64[0 1 3 4 5 0;
                                                  1 0 3 4 5 1;
                                                  3 3 0 0 0 1;
                                                  4 4 0 0 0 1;
                                                  5 5 0 0 0 1;
                                                  0 1 1 1 1 0]), 1)
            cb_set = Set{DiameterSimplex{2, Float64}}()
            for sx in coboundary(st, DiameterSimplex{2}(st, 1.0, (2, 1), 1), 1)
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
            st = ReductionState{3}(sparse(dist), 5)

            for dim in 1:5
                cob = DiameterSimplex{3, Float64}[]
                for sx in coboundary(st, DiameterSimplex{3}(1.0, 10, 1), dim)
                    push!(cob, sx)
                end
                # should be reverse sorted?
                @test issorted(index.(cob), rev=true)
                @test length(cob) == 100 - dim - 1
            end
        end

        @testset "diameters" begin
            dist = torus(25)
            st = ReductionState{2}(dist, 3)
            for dim in 1:3
                diameter = diam(st, dim+1:-1:1)
                sx = DiameterSimplex{2}(st, diameter, 1, 1)
                for coface in coboundary(st, sx, dim)
                    @test diam(coface) == diam(st, vertices(st, coface, dim+1))
                end
            end
        end
    end

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
            @test isnothing(pop_pivot!(col))
        end
    end

    @testset "compute_0_dim_pairs!" begin
        @testset "dense Int" begin
            dist = [0 1 2;
                    1 0 3;
                    2 3 0]
            state = ReductionState{2}(dist, 0)
            simplices = DiameterSimplex{2, Int}[]
            critical_edges = DiameterSimplex{2, Int}[]
            res = compute_0_dim_pairs!(state, simplices, critical_edges)

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
            state = ReductionState{3}(dist, 0)
            simplices = DiameterSimplex{3, Float64}[]
            critical_edges = DiameterSimplex{3, Float64}[]
            res = compute_0_dim_pairs!(state, simplices, critical_edges)

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
        state = ReductionState{2}(dist, 1)
        simplices = DiameterSimplex{2, Float64}[]
        columns = DiameterSimplex{2, Float64}[]
        compute_0_dim_pairs!(state, simplices, columns)
        matrices = ReductionMatrices(state, 1)

        res = compute_pairs!(matrices, columns)

        @test res == [(1.0, 2.0)]
    end

    @testset "ripserer" begin
        @testset "icosahedron" begin
            icosahedron = Float64[0  1  2  2  1  2  1  1  2  2  1  3;
                                  1  0  3  2  1  1  2  1  2  1  2  2;
                                  2  3  0  1  2  2  1  2  1  2  1  1;
                                  2  2  1  0  3  2  1  1  2  1  2  1;
                                  1  1  2  3  0  1  2  2  1  2  1  2;
                                  2  1  2  2  1  0  3  2  1  1  2  1;
                                  1  2  1  1  2  3  0  1  2  2  1  2;
                                  1  1  2  1  2  2  1  0  3  1  2  2;
                                  2  2  1  2  1  1  2  3  0  2  1  1;
                                  2  1  2  1  2  1  2  1  2  0  3  1;
                                  1  2  1  2  1  2  1  2  1  3  0  2;
                                  3  2  1  1  2  1  2  2  1  1  2  0]

            res = ripserer(icosahedron, 2)
            @test res[1] == [fill((0.0, 1.0), 11); (0.0, Inf)]
            @test isempty(res[2])
            @test res[3] == [(1.0, 2.0)]
        end
        @testset "torus 16" begin
            dim0, dim1, dim2 = ripserer(torus(16), 2)

            @test length(dim0) == 16

            @test all(x -> first(x) ≈ 0.5, dim1)
            @test sum(x -> last(x) ≈ 1, dim1) == 2
            @test sum(x -> isapprox(last(x), 0.71, atol=0.1), dim1) == 15

            @test last(only(dim2)) == 1
        end
        @testset "torus 100" begin
            dim0, dim1 = ripserer(torus(100), 1)

            @test length(dim0) == 100

            deaths = sort(last.(dim1))
            @test deaths[end] ≈ 0.8
            @test deaths[end-1] ≈ 0.8
            @test deaths[end-2] < 0.5
        end
        @testset "cycle" begin
            cycle = [0 1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2 1;
                     1 0 1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2;
                     2 1 0 1 2 3 4 5 6 7 8 9 8 7 6 5 4 3;
                     3 2 1 0 1 2 3 4 5 6 7 8 9 8 7 6 5 4;
                     4 3 2 1 0 1 2 3 4 5 6 7 8 9 8 7 6 5;
                     5 4 3 2 1 0 1 2 3 4 5 6 7 8 9 8 7 6;
                     6 5 4 3 2 1 0 1 2 3 4 5 6 7 8 9 8 7;
                     7 6 5 4 3 2 1 0 1 2 3 4 5 6 7 8 9 8;
                     8 7 6 5 4 3 2 1 0 1 2 3 4 5 6 7 8 9;
                     9 8 7 6 5 4 3 2 1 0 1 2 3 4 5 6 7 8;
                     8 9 8 7 6 5 4 3 2 1 0 1 2 3 4 5 6 7;
                     7 8 9 8 7 6 5 4 3 2 1 0 1 2 3 4 5 6;
                     6 7 8 9 8 7 6 5 4 3 2 1 0 1 2 3 4 5;
                     5 6 7 8 9 8 7 6 5 4 3 2 1 0 1 2 3 4;
                     4 5 6 7 8 9 8 7 6 5 4 3 2 1 0 1 2 3;
                     3 4 5 6 7 8 9 8 7 6 5 4 3 2 1 0 1 2;
                     2 3 4 5 6 7 8 9 8 7 6 5 4 3 2 1 0 1;
                     1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2 1 0]

            dim0, dim1, dim2, dim3, dim4 = ripserer(cycle, 4)
            @test dim0 == [fill((0, 1), size(cycle, 1) - 1); (0, typemax(Int))]
            @test dim1 == [(1, 6)]
            @test dim2 == fill((6, 7), 5)
            @test dim3 == [(7, 8)]
            @test dim4 == []
        end
    end
end
