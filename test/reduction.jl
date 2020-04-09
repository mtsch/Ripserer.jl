using Ripserer:
    CompressedSparseMatrix, add_column!,
    Column, coboundary, pop_pivot!,
    compute_0_dim_pairs!,
    ReductionMatrix, compute_pairs!

@testset "reduction" begin
    @testset "CompressedSparseMatrix" begin
        csm = CompressedSparseMatrix{Int}()
        @test length(csm) == 0

        add_column!(csm)
        push!(csm, 1)
        push!(csm, 2)
        push!(csm, 3)
        push!(csm, 4)

        add_column!(csm)
        push!(csm, 0)
        push!(csm, 0)
        push!(csm, 0)

        add_column!(csm)

        add_column!(csm)
        push!(csm, 1)

        @test length(csm[1]) == 4
        @test collect(csm[1]) == [1, 2, 3, 4]

        @test length(csm[2]) == 3
        @test all(iszero, csm[2])

        @test length(csm[3]) == 0
        @test eltype(collect(csm[3])) === Int

        @test length(csm[4]) == 1
        @test first(csm[4]) == 1

        @test length(csm) == 4
    end

    @testset "a Column is a heap." begin
        column = Column{3, Float64}()
        @test isempty(column)

        push!(column, Simplex{3}(3.0, 1, 1))
        push!(column, Simplex{3}(3.0, 2, 1))
        push!(column, Simplex{3}(1.0, 3, 1))
        push!(column, Simplex{3}(2.0, 4, 1))
        push!(column, Simplex{3}(2.0, 4, 1))
        @test length(column) == 5
        @test !isempty(column)

        @test top(column) == Simplex{3}(1.0, 3, 1)
        @test pop!(column) == Simplex{3}(1.0, 3, 1)
        @test length(column) == 4
        @test pop!(column) == Simplex{3}(2.0, 4, 1)
        @test pop!(column) == Simplex{3}(2.0, 4, 1)
        @test pop!(column) == Simplex{3}(3.0, 2, 1)
        @test pop!(column) == Simplex{3}(3.0, 1, 1)
        @test isempty(column)

        @test_throws MethodError push!(column, Simplex{2}(3.0, 1, 1))
        @test_throws MethodError push!(column, Simplex{3}(3, 1, 1))
    end

    @testset "pop_pivot!" begin
        @testset "single element" begin
            col = Column{2, Float64}()
            push!(col, Simplex{2}(2.0, 1, 1))
            push!(col, Simplex{2}(2.0, 1, 1))
            push!(col, Simplex{2}(2.0, 1, 1))
            push!(col, Simplex{2}(2.0, 1, 1))
            push!(col, Simplex{2}(2.0, 1, 1))

            @test pop_pivot!(col) == Simplex{2}(2.0, 1, 1)
            @test isempty(col)

            col = Column{3, Float64}()
            push!(col, Simplex{3}(2.0, 1, 1))
            push!(col, Simplex{3}(2.0, 1, 1))
            push!(col, Simplex{3}(2.0, 1, 1))

            @test isnothing(pop_pivot!(col))
            @test isempty(col)
        end

        @testset "multiple" begin
            col = Column{5, Float64}()
            push!(col, Simplex{5}(1.0, 2, 3))
            push!(col, Simplex{5}(2.0, 3, 4))
            push!(col, Simplex{5}(1.0, 2, 2))
            push!(col, Simplex{5}(3.0, 1, 2))
            push!(col, Simplex{5}(2.0, 3, 1))
            push!(col, Simplex{5}(4.0, 4, 4))
            push!(col, Simplex{5}(4.0, 4, 4))
            push!(col, Simplex{5}(4.0, 4, 4))

            @test pop_pivot!(col) == Simplex{5}(3.0, 1, 2)
            @test pop_pivot!(col) == Simplex{5}(4.0, 4, 2)
            @test isnothing(pop_pivot!(col))
            @test isnothing(pop_pivot!(col))
        end
    end

    @testset "compute_0_dim_pairs!" begin
        @testset "dense Int" begin
            dist = [0 1 2;
                    1 0 3;
                    2 3 0]
            scx = RipsComplex{2}(dist, 0)
            critical_edges = Simplex{2, Int}[]
            res = compute_0_dim_pairs!(scx, critical_edges)

            @test res == [(0, 1),
                          (0, 2),
                          (0, typemax(Int))]
            @test critical_edges == [Simplex{2}(3, 3, 1)]
        end

        #=
        @testset "sparse Float64" begin
            dist = sparse(Float64[0 2 0 0 5 0;
                                  2 0 4 6 0 0;
                                  0 4 0 3 0 0;
                                  0 6 3 0 1 0;
                                  5 0 0 1 0 0;
                                  0 0 0 0 0 0])
            scx = RipsComplex{3}(dist, 0)
            simplices = Simplex{3, Float64}[]
            critical_edges = Simplex{3, Float64}[]
            res = compute_0_dim_pairs!(scx, simplices, critical_edges)

            @test res == [(0.0, 1.0),
                          (0.0, 2.0),
                          (0.0, 3.0),
                          (0.0, 4.0),
                          (0.0, Inf),
                          (0.0, Inf)]
            @test critical_edges == [Simplex{3}(6.0, 5, 1),
                                     Simplex{3}(5.0, 7, 1)]
        end
        =#
    end

    @testset "ripserer" begin
        @testset "full matrix, no threshold" begin
            @testset "icosahedron" begin
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

                dim0_7, dim1_7, dim2_7, dim3_7, dim4_7 = ripserer(cycle, 4, 7)
                @test all(dim0 .== dim0_7)
                @test all(dim1 .== dim1_7)
                @test all(dim2 .== dim2_7)
                @test all(dim3 .== dim3_7)
                @test all(dim4 .== dim4_7)
            end

            @testset "projective plane (modulus)" begin
                # taken from ripser/examples
                projective_plane = [0 1 1 1 1 1 1 1 1 2 2 2 2;
                                    1 0 2 2 2 1 2 1 2 1 2 2 2;
                                    1 2 0 2 2 2 1 2 1 1 2 2 2;
                                    1 2 2 0 2 1 2 2 1 2 2 2 1;
                                    1 2 2 2 0 2 1 1 2 2 2 2 1;
                                    1 1 2 1 2 0 2 2 2 1 1 2 1;
                                    1 2 1 2 1 2 0 2 2 1 1 2 1;
                                    1 1 2 2 1 2 2 0 2 1 2 1 1;
                                    1 2 1 1 2 2 2 2 0 1 2 1 1;
                                    2 1 1 2 2 1 1 1 1 0 1 1 2;
                                    2 2 2 2 2 1 1 2 2 1 0 2 1;
                                    2 2 2 2 2 2 2 1 1 1 2 0 1;
                                    2 2 2 1 1 1 1 1 1 2 1 1 0]

                _, dim1_2, dim2_2 = ripserer(projective_plane, 2)
                _, dim1_3, dim2_3 = ripserer(projective_plane, 2, 3)
                @test dim1_2 == [(1, 2)]
                @test dim2_2 == [(1, 2)]
                @test isempty(dim1_3)
                @test isempty(dim2_3)
            end
        end
    end
end
