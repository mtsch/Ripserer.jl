using Ripserer: Infinity, ∞, edges, is_distance_matrix, distances

@testset "convenience" begin
    @testset "edges dense" begin
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
    @testset "edges sparse" begin
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
    @testset "n edges dense" begin
        dist = rand_dist_matrix(100)
        n_edges = binomial(size(dist, 1), 2)
        @test length(edges(dist)) == n_edges
    end
    @testset "n edges sparse" begin
        for _ in 1:10
            dist = rand_dist_matrix(100, 0.5)
            n_edges = nnz(dist) ÷ 2
            @test length(edges(dist)) == n_edges
        end
    end
    @testset "is_distance_matrix" begin
        for f in (identity, sparse)
            @test is_distance_matrix(f([0 1 0;
                                        1 0 2;
                                        0 2 0]))
            @test is_distance_matrix(f([0 1 0 0;
                                        1 0 2 0;
                                        0 2 0 0;
                                        0 0 0 0]))
            @test !is_distance_matrix(f([1 1 0;
                                         1 0 2;
                                         0 2 0]))
            @test !is_distance_matrix(f([0 2 0;
                                         1 0 2;
                                         0 2 0]))
            @test !is_distance_matrix(f([0 1 0 0;
                                         1 0 2 0;
                                         0 2 0 0]))
        end
    end
    @testset "distances" begin
        points = [(0, 0), (0, 1), (1, 1), (1, 0)]
        @test distances(Euclidean(), points) ≈ [ 0  1 √2  1;
                                                 1  0  1 √2;
                                                √2  1  0  1;
                                                 1 √2  1  0]

        @test distances(Cityblock(), points) == [0 1 2 1;
                                                 1 0 1 2;
                                                 2 1 0 1;
                                                 1 2 1 0]
    end
    @testset "infinty" begin
        @test ∞ > 0.0
        @test typemax(Int) < ∞
        @test ∞ == Inf
        @test Inf == ∞
        @test !(∞ > Inf)
        @test !(∞ < Inf)
        @test !(∞ > NaN)
        @test !(∞ < NaN)
        @test !(∞ > ∞)
        @test !(∞ < ∞)
        @test ∞ > "infinity plus one"
        @test "infinity plus two" < ∞
        @test ∞ > ["infinity", +, 3]
        @test ismissing(∞ > missing)
        @test ismissing(∞ < missing)
        @test !isless(∞, missing)
        @test isless(missing, ∞)
        @test ∞ ≈ Inf
        @test 1 ≉ Inf

        @test sprint(print, ∞) == "∞"
    end
end
