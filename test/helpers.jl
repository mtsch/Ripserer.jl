using Ripserer:
    edges, is_distance_matrix, apply_threshold,
    isprime, Binomial, ReductionState, n_vertices, dist, is_connected,
    CompressedSparseMatrix, add_column!

@testset "helpers" begin
    @testset "distancematrix" begin
        @testset "edges" begin
            @testset "dense" begin
                dist = [0 1 2 3;
                        1 0 4 5;
                        2 4 0 4;
                        3 5 4 0]
                res = edges(dist)
                @test res[1] == (1, (2, 1))
                @test res[2] == (2, (3, 1))
                @test res[3] == (3, (4, 1))
                @test res[4] == (4, (4, 3))
                @test res[5] == (4, (3, 2))
                @test res[6] == (5, (4, 2))
            end

            @testset "sparse" begin
                dist = sparse([0 1 0 3;
                               1 0 4 5;
                               0 4 0 4;
                               3 5 4 0])
                res = edges(dist)
                @test res[1] == (1, (2, 1))
                @test res[2] == (3, (4, 1))
                @test res[3] == (4, (4, 3))
                @test res[4] == (4, (3, 2))
                @test res[5] == (5, (4, 2))
            end

            @testset "number of edges" begin
                for i in 1:10
                    dist = rand_dist_matrix(100, 0.5)
                    n_edges = nnz(dist) ÷ 2
                    @test length(edges(dist)) == n_edges
                end
            end
        end

        @testset "is_distance_matrix" begin
            @test is_distance_matrix([0 1 0;
                                      1 0 2;
                                      0 2 0])
            @test is_distance_matrix([0 1 0 0;
                                      1 0 2 0;
                                      0 2 0 0;
                                      0 0 0 0])
            @test !is_distance_matrix([1 1 0;
                                       1 0 2;
                                       0 2 0])
            @test !is_distance_matrix([0 2 0;
                                       1 0 2;
                                       0 2 0])
            @test !is_distance_matrix([0 1 0 0;
                                       1 0 2 0;
                                       0 2 0 0])
        end

        @testset "apply_threshold" begin
            for t in 0:0.2:1
                dist_d = apply_threshold(rand_dist_matrix(10), t)
                dist_s = apply_threshold(rand_dist_matrix(10, 0.5), t)

                @test maximum(dist_d) ≤ t
                @test maximum(dist_s) ≤ t
                @test is_distance_matrix(dist_d)
                @test is_distance_matrix(dist_s)
                @test issparse(dist_d)
                @test issparse(dist_s)
            end
        end
    end

    @testset "state" begin
        @testset "isprime" begin
            @test !isprime(1)
            @test isprime(2)
            @test isprime(3)
            @test !isprime(4)
            @test isprime(5)
            @test !isprime(6)
            @test isprime(7)
            @test !isprime(8)
            @test !isprime(9)
            @test !isprime(10)
            @test isprime(11)
        end

        @testset "Binomial" begin
            bin = Binomial(10, 15)
            for n in 1:10, k in 1:15
                @test bin(n, k) == binomial(n, k)
            end
        end

        @testset "ReductionState" begin
            st = ReductionState{3}([0 1 2 0;
                                    1 0 3 0;
                                    2 3 0 0;
                                    0 0 0 0], 1)
            @test n_vertices(st) == 4
            @test dist(st, 3, 3) == 0
            @test dist(st, 1, 2) == 1
            @test dist(st, 1, 2) == 1
            @test is_connected(st, 1, 2)
            @test !is_connected(st, 3, 3)
            @test !is_connected(st, 2, 4)
        end
    end

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
end
