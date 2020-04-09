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
                res = ripserer(icosahedron, dim_max=2)
                @test res[1] == [fill((0.0, 1.0), 11); (0.0, Inf)]
                @test isempty(res[2])
                @test res[3] == [(1.0, 2.0)]
            end

            @testset "torus 16" begin
                d0, d1, d2 = ripserer(torus(16), dim_max=2)

                @test length(d0) == 16

                @test all(x -> first(x) ≈ 0.5, d1)
                @test sum(x -> last(x) ≈ 1, d1) == 2
                @test sum(x -> isapprox(last(x), 0.71, atol=0.1), d1) == 15

                @test last(only(d2)) == 1
            end

            @testset "torus 100" begin
                d0, d1 = ripserer(torus(100), dim_max=1)

                @test length(d0) == 100

                deaths = sort(last.(d1))
                @test deaths[end] ≈ 0.8
                @test deaths[end-1] ≈ 0.8
                @test deaths[end-2] < 0.5
            end

            @testset "cycle" begin
                d0, d1, d2, d3, d4 = ripserer(cycle, dim_max=4)
                @test d0 == [fill((0, 1), size(cycle, 1) - 1); (0, typemax(Int))]
                @test d1 == [(1, 6)]
                @test d2 == fill((6, 7), 5)
                @test d3 == [(7, 8)]
                @test d4 == []

                d0_7, d1_7, d2_7, d3_7, d4_7 = ripserer(cycle, dim_max=4, modulus=7)
                @test all(d0 .== d0_7)
                @test all(d1 .== d1_7)
                @test all(d2 .== d2_7)
                @test all(d3 .== d3_7)
                @test all(d4 .== d4_7)
            end

            @testset "projective plane (modulus)" begin
                _, d1_2, d2_2 = ripserer(projective_plane, dim_max=2)
                _, d1_3, d2_3 = ripserer(projective_plane, dim_max=2, modulus=3)
                @test d1_2 == [(1, 2)]
                @test d2_2 == [(1, 2)]
                @test isempty(d1_3)
                @test isempty(d2_3)
            end
        end

        @testset "full matrix, with threshold" begin
            @testset "icosahedron, high thresh" begin
                res = ripserer(icosahedron, threshold=2, dim_max=2)
                @test res[1] == [fill((0.0, 1.0), 11); (0.0, Inf)]
                @test isempty(res[2])
                @test res[3] == [(1.0, 2.0)]
            end

            @testset "icosahedron, med thresh" begin
                res = ripserer(icosahedron, dim_max=2, threshold=1)
                @test res[1] == [fill((0.0, 1.0), 11); (0.0, Inf)]
                @test isempty(res[2])
                @test res[3] == [(1.0, Inf)]
            end

            @testset "icosahedron, low thresh" begin
                res = ripserer(icosahedron, dim_max=2, threshold=0.5)
                @test res[1] == fill((0.0, Inf), 12)
                @test isempty(res[2])
                @test isempty(res[3])
            end

            @testset "torus 16, high threshold" begin
                d0, d1, d2 = ripserer(torus(16), dim_max=2, threshold=2)

                @test length(d0) == 16

                @test all(x -> first(x) ≈ 0.5, d1)
                @test sum(x -> last(x) ≈ 1, d1) == 2
                @test sum(x -> isapprox(last(x), 0.71, atol=0.1), d1) == 15

                @test last(only(d2)) == 1
            end

            @testset "torus 16, med threshold" begin
                d0, d1, d2 = ripserer(torus(16), dim_max=2, threshold=0.9)

                @test length(d0) == 16

                @test all(x -> first(x) ≈ 0.5, d1)
                @test sum(x -> last(x) == Inf, d1) == 2
                @test sum(x -> isapprox(last(x), 0.71, atol=0.1), d1) == 15

                @test last(only(d2)) == Inf
            end

            @testset "torus 16, low threshold" begin
                d0, d1, d2 = ripserer(torus(16), dim_max=2, threshold=0.5)

                @test length(d0) == 16

                @test all(x -> first(x) ≈ 0.5, d1)
                @test all(x -> last(x) == Inf, d1)

                @test isempty(d2)
            end

            @testset "projective plane (modulus), med threshold" begin
                _, d1_2, d2_2 = ripserer(projective_plane,
                                         dim_max=2, threshold=1)
                _, d1_3, d2_3 = ripserer(projective_plane,
                                         dim_max=2, modulus=3, threshold=1)
                @test d1_2 == [(1, typemax(Int))]
                @test d2_2 == [(1, typemax(Int))]
                @test isempty(d1_3)
                @test isempty(d2_3)
            end
        end
    end
end
