using Ripserer: DiameterSimplex, copy_edges!, edges, is_distance_matrix, apply_threshold
using SparseArrays

function rand_dist_matrix(n, sparse::Bool=false)
    A = rand(n, n)
    A .+= A'
    A -= Diagonal(A)
    A
end
function rand_dist_matrix(n, sparse)
    A = sprand(n, n, sparse/2)
    A .+= A'
    A -= Diagonal(A)
    A
end

@testset "distancematrix" begin
    @testset "copy_edges!: dense" begin
        buff = DiameterSimplex{2, Float64}[]
        dist = [0 1 2 3;
                1 0 4 5;
                2 4 0 4;
                3 5 4 0]
        copy_edges!(buff, dist, binomial)
        @test buff[1] == DiameterSimplex{2}(1.0, 1, 1)
        @test buff[2] == DiameterSimplex{2}(2.0, 2, 1)
        @test buff[3] == DiameterSimplex{2}(3.0, 4, 1)
        @test buff[4] == DiameterSimplex{2}(4.0, 6, 1)
        @test buff[5] == DiameterSimplex{2}(4.0, 3, 1)
        @test buff[6] == DiameterSimplex{2}(5.0, 5, 1)
    end

    @testset "copy_edges!: sparse" begin
        buff = DiameterSimplex{3, Int}[]
        dist = sparse([0 1 0 3;
                       1 0 4 5;
                       0 4 0 4;
                       3 5 4 0])
        copy_edges!(buff, dist, binomial)
        @test buff[1] == DiameterSimplex{3}(1, 1, 1)
        @test buff[2] == DiameterSimplex{3}(3, 4, 1)
        @test buff[3] == DiameterSimplex{3}(4, 6, 1)
        @test buff[4] == DiameterSimplex{3}(4, 3, 1)
        @test buff[5] == DiameterSimplex{3}(5, 5, 1)
    end

    @testset "copy_edges!: number of edges" begin
        for i in 1:10
            dist = rand_dist_matrix(100, 0.5)
            n_edges = nnz(dist) ÷ 2
            buff = DiameterSimplex{7, Float64}[]
            @test length(copy_edges!(buff, dist, binomial)) == n_edges
        end
    end

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
