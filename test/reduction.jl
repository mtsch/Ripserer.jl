using Ripserer:
    CompressedSparseMatrix, add_column!,
    Column, pop_pivot!,
    zeroth_intervals,
    ReductionMatrix

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

    @testset "pop_pivot!" begin
        @testset "single element" begin
            col = Column{Simplex{1, 2, Float64}}()
            push!(col, Simplex{1, 2}(2.0, 1, 1))
            push!(col, Simplex{1, 2}(2.0, 1, 1))
            push!(col, Simplex{1, 2}(2.0, 1, 1))
            push!(col, Simplex{1, 2}(2.0, 1, 1))
            push!(col, Simplex{1, 2}(2.0, 1, 1))

            @test pop_pivot!(col) == Simplex{1, 2}(2.0, 1, 1)
            @test isempty(col)

            col = Column{Simplex{2, 3, Float64}}()
            push!(col, Simplex{2, 3}(2.0, 1, 1))
            push!(col, Simplex{2, 3}(2.0, 1, 1))
            push!(col, Simplex{2, 3}(2.0, 1, 1))

            @test isnothing(pop_pivot!(col))
            @test isempty(col)
        end
        @testset "multiple" begin
            col = Column{Simplex{3, 5, Float64}}()
            push!(col, Simplex{3, 5}(1.0, 2, 3))
            push!(col, Simplex{3, 5}(2.0, 3, 4))
            push!(col, Simplex{3, 5}(1.0, 2, 2))
            push!(col, Simplex{3, 5}(3.0, 1, 2))
            push!(col, Simplex{3, 5}(2.0, 3, 1))
            push!(col, Simplex{3, 5}(4.0, 4, 4))
            push!(col, Simplex{3, 5}(4.0, 4, 4))
            push!(col, Simplex{3, 5}(4.0, 4, 4))
            push!(col, Simplex{3, 5}(5.0, 4, 4))
            push!(col, Simplex{3, 5}(5.0, 4, 1))

            @test pop_pivot!(col) == Simplex{3, 5}(3.0, 1, 2)
            @test pop_pivot!(col) == Simplex{3, 5}(4.0, 4, 2)
            @test isnothing(pop_pivot!(col))
            @test isnothing(pop_pivot!(col))
        end
    end

    @testset "compute_0_dim_pairs!" begin
        @testset "dense" begin
            dist = [0 1 2;
                    1 0 3;
                    2 3 0]
            flt = RipsFiltration(dist, threshold=3)
            res, columns, simplices = zeroth_intervals(flt)

            @test isnothing(simplices)
            @test res == [(0, 1),
                          (0, 2),
                          (0, ∞)]
            @test columns == [Simplex{1, 2}(3, 3, 1)]
        end
        @testset "sparse" begin
            dist = [0 1 2;
                    1 0 3;
                    2 3 0]
            flt = SparseRipsFiltration(dist)
            res, columns, simplices = zeroth_intervals(flt)

            @test simplices == [Simplex{1, 2}(1, 1, 1),
                                Simplex{1, 2}(2, 2, 1)]
            @test res == [(0, 1),
                          (0, 2),
                          (0, ∞)]
            @test isempty(columns)
        end
    end

    @testset "ripserer" begin
        @testset "full matrix, no threshold" begin
            @testset "icosahedron" begin
                res = ripserer(icosahedron, dim_max=2)
                @test res[1] == [fill((0.0, 1.0), 11); (0.0, ∞)]
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
                @test d0 == [fill((0, 1), size(cycle, 1) - 1); (0, ∞)]
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
                @test res[1] == [fill((0.0, 1.0), 11); (0.0, ∞)]
                @test isempty(res[2])
                @test res[3] == [(1.0, 2.0)]
            end
            @testset "icosahedron, med thresh" begin
                res = ripserer(icosahedron, dim_max=2, threshold=1)
                @test res[1] == [fill((0.0, 1.0), 11); (0.0, ∞)]
                @test isempty(res[2])
                @test res[3] == [(1.0, ∞)]
            end
            @testset "icosahedron, low thresh" begin
                res = ripserer(icosahedron, dim_max=2, threshold=0.5)
                @test res[1] == fill((0.0, ∞), 12)
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
                @test sum(x -> last(x) == ∞, d1) == 2
                @test sum(x -> isapprox(last(x), 0.71, atol=0.1), d1) == 15

                @test last(only(d2)) == ∞
            end
            @testset "torus 16, low threshold" begin
                d0, d1, d2 = ripserer(torus(16), dim_max=2, threshold=0.5)

                @test length(d0) == 16

                @test all(x -> first(x) ≈ 0.5, d1)
                @test all(x -> last(x) == ∞, d1)

                @test isempty(d2)
            end
            @testset "projective plane (modulus), med threshold" begin
                _, d1_2, d2_2 = ripserer(projective_plane,
                                         dim_max=2, threshold=1)
                _, d1_3, d2_3 = ripserer(projective_plane,
                                         dim_max=2, modulus=3, threshold=1)
                @test d1_2 == [(1, ∞)]
                @test d2_2 == [(1, ∞)]
                @test isempty(d1_3)
                @test isempty(d2_3)
            end
        end

        @testset "sparse matrix" begin
            @testset "icosahedron" begin
                flt = SparseRipsFiltration(icosahedron, threshold=2)
                res = ripserer(flt, dim_max=2)
                @test res[1] == [fill((0.0, 1.0), 11); (0.0, ∞)]
                @test isempty(res[2])
                @test res[3] == [(1.0, 2.0)]
            end
            @testset "torus 16" begin
                dists = sparse(torus(16))
                SparseArrays.fkeep!(dists, (_, _, v) -> v ≤ 1)

                d0, d1, d2 = ripserer(dists, dim_max=2)

                @test length(d0) == 16

                @test all(x -> first(x) ≈ 0.5, d1)
                @test sum(x -> last(x) ≈ 1, d1) == 2
                @test sum(x -> isapprox(last(x), 0.71, atol=0.1), d1) == 15

                @test last(only(d2)) == 1
            end
            @testset "projective plane (modulus), med threshold" begin
                dists = sparse(projective_plane)
                SparseArrays.fkeep!(dists, (_, _, v) -> v ≤ 2)

                _, d1_2, d2_2 = ripserer(dists, dim_max=2, threshold=1)
                _, d1_3, d2_3 = ripserer(dists, dim_max=2, modulus=3, threshold=1)
                @test d1_2 == [(1, ∞)]
                @test d2_2 == [(1, ∞)]
                @test isempty(d1_3)
                @test isempty(d2_3)
            end
        end
    end
end
