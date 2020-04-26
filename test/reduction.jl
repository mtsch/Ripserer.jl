using Ripserer:
    ReductionMatrix, insert_column!, has_column,
    Column, pop_pivot!,
    zeroth_intervals,
    ReductionMatrix

@testset "reduction" begin
    @testset "ReductionMatrix" begin
        rm = ReductionMatrix{Int}()
        @test length(rm) == 0

        insert_column!(rm, 3)
        push!(rm, 1)
        push!(rm, 2)
        push!(rm, 3)
        push!(rm, 4)

        insert_column!(rm, 10)
        push!(rm, 0)
        push!(rm, 0)
        push!(rm, 0)

        insert_column!(rm, 1)

        insert_column!(rm, 15)
        push!(rm, 1)

        @test has_column(rm, 3)
        @test collect(rm[3]) == [1, 2, 3, 4]

        @test has_column(rm, 10)
        @test length(rm[10]) == 3
        @test all(iszero, rm[10])

        @test has_column(rm, 1)
        @test length(rm[1]) == 0
        @test eltype(collect(rm[1])) === Int

        @test has_column(rm, 15)
        @test length(rm[15]) == 1
        @test first(rm[15]) == 1

        @test !has_column(rm, 2)
        @test !has_column(rm, 100)
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

            @test !isnothing(simplices)
            @test res == [PersistenceInterval(0, 1),
                          PersistenceInterval(0, 2),
                          PersistenceInterval(0, ∞)]
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
            @test res == [PersistenceInterval(0, 1),
                          PersistenceInterval(0, 2),
                          PersistenceInterval(0, ∞)]
            @test isempty(columns)
        end
    end

    @testset "ripserer" begin
        @testset "full matrix, no threshold" begin
            @testset "icosahedron" begin
                res = ripserer(icosahedron, dim_max=2)
                @test res[1] == [fill(PersistenceInterval(0.0, 1.0), 11);
                                 PersistenceInterval(0.0, ∞)]
                @test isempty(res[2])
                @test res[3] == [PersistenceInterval(1.0, 2.0)]
            end
            @testset "torus 16" begin
                d0, d1, d2 = ripserer(torus(16), dim_max=2)

                @test length(d0) == 16

                @test all(x -> birth(x) ≈ 0.5, d1)
                @test count(x -> death(x) ≈ 1, d1) == 2
                @test count(x -> isapprox(death(x), 0.71, atol=0.1), d1) == 15

                @test death(only(d2)) == 1
            end
            @testset "torus 100" begin
                d0, d1 = ripserer(torus(100), dim_max=1)

                @test length(d0) == 100

                deaths = sort(death.(d1))
                @test deaths[end] ≈ 0.8
                @test deaths[end-1] ≈ 0.8
                @test deaths[end-2] < 0.5
            end
            @testset "cycle" begin
                d0, d1, d2, d3, d4 = ripserer(cycle, dim_max=4)
                @test d0 == [fill(PersistenceInterval(0, 1), size(cycle, 1) - 1);
                             PersistenceInterval(0, ∞)]
                @test d1 == [PersistenceInterval(1, 6)]
                @test d2 == fill(PersistenceInterval(6, 7), 5)
                @test d3 == [PersistenceInterval(7, 8)]
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
                @test d1_2 == [PersistenceInterval(1, 2)]
                @test d2_2 == [PersistenceInterval(1, 2)]
                @test isempty(d1_3)
                @test isempty(d2_3)
            end
        end

        @testset "full matrix, with threshold" begin
            @testset "icosahedron, high threshold" begin
                res = ripserer(icosahedron, threshold=2, dim_max=2)
                @test res[1] == [fill(PersistenceInterval(0.0, 1.0), 11);
                                 PersistenceInterval(0.0, ∞)]
                @test isempty(res[2])
                @test res[3] == [PersistenceInterval(1.0, 2.0)]
            end
            @testset "icosahedron, med threshold" begin
                res = ripserer(icosahedron, dim_max=2, threshold=1)
                @test res[1] == [fill(PersistenceInterval(0.0, 1.0), 11);
                                 PersistenceInterval(0.0, ∞)]
                @test isempty(res[2])
                @test res[3] == [PersistenceInterval(1.0, ∞)]
            end
            @testset "icosahedron, low threshold" begin
                res = ripserer(icosahedron, dim_max=2, threshold=0.5)
                @test res[1] == fill(PersistenceInterval(0.0, ∞), 12)
                @test isempty(res[2])
                @test isempty(res[3])
            end
            @testset "torus 16, high threshold" begin
                d0, d1, d2 = ripserer(torus(16), dim_max=2, threshold=2)

                @test length(d0) == 16

                @test all(x -> birth(x) ≈ 0.5, d1)
                @test count(x -> death(x) ≈ 1, d1) == 2
                @test count(x -> isapprox(death(x), 0.71, atol=0.1), d1) == 15

                @test death(only(d2)) == 1
            end
            @testset "torus 16, med threshold" begin
                d0, d1, d2 = ripserer(torus(16), dim_max=2, threshold=0.9)

                @test length(d0) == 16

                @test all(x -> birth(x) ≈ 0.5, d1)
                @test count(x -> death(x) == ∞, d1) == 2
                @test count(x -> isapprox(death(x), 0.71, atol=0.1), d1) == 15

                @test last(only(d2)) == ∞
            end
            @testset "torus 16, low threshold" begin
                d0, d1, d2 = ripserer(torus(16), dim_max=2, threshold=0.5)

                @test length(d0) == 16

                @test all(x -> birth(x) ≈ 0.5, d1)
                @test all(x -> death(x) == ∞, d1)

                @test isempty(d2)
            end
            @testset "projective plane (modulus), med threshold" begin
                _, d1_2, d2_2 = ripserer(projective_plane,
                                         dim_max=2, threshold=1)
                _, d1_3, d2_3 = ripserer(projective_plane,
                                         dim_max=2, modulus=3, threshold=1)
                @test d1_2 == [PersistenceInterval(1, ∞)]
                @test d2_2 == [PersistenceInterval(1, ∞)]
                @test isempty(d1_3)
                @test isempty(d2_3)
            end
        end

        @testset "sparse matrix" begin
            @testset "icosahedron" begin
                flt = SparseRipsFiltration(icosahedron, threshold=2)
                res = ripserer(flt, dim_max=2)
                @test res[1] == [fill(PersistenceInterval(0.0, 1.0), 11);
                                 PersistenceInterval(0.0, ∞)]
                @test isempty(res[2])
                @test res[3] == [PersistenceInterval(1.0, 2.0)]
            end
            @testset "torus 16" begin
                dists = sparse(torus(16))
                SparseArrays.fkeep!(dists, (_, _, v) -> v ≤ 1)

                d0, d1, d2 = ripserer(dists, dim_max=2)

                @test length(d0) == 16

                @test all(x -> birth(x) ≈ 0.5, d1)
                @test count(x -> death(x) ≈ 1, d1) == 2
                @test count(x -> isapprox(death(x), 0.71, atol=0.1), d1) == 15

                @test last(only(d2)) == 1
            end
            @testset "projective plane (modulus), med threshold" begin
                dists = sparse(projective_plane)
                SparseArrays.fkeep!(dists, (_, _, v) -> v ≤ 2)

                _, d1_2, d2_2 = ripserer(dists, dim_max=2, threshold=1)
                _, d1_3, d2_3 = ripserer(dists, dim_max=2, modulus=3, threshold=1)
                @test d1_2 == [PersistenceInterval(1, ∞)]
                @test d2_2 == [PersistenceInterval(1, ∞)]
                @test isempty(d1_3)
                @test isempty(d2_3)
            end
        end
    end
end
