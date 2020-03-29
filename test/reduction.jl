using Ripserer: compute_0_dim_pairs!, copy_edges!

@testset "compute_0_dim_pairs!" begin
    @testset "dense Int" begin
        dist = [0 1 2;
                1 0 3;
                2 3 0]
        edges = copy_edges!(DiameterSimplex{2, Float64}[], dist, binomial)
        critical_edges = DiameterSimplex{2, Int64}[]

        res = compute_0_dim_pairs!(critical_edges, dist, binomial)

        @test res == [(0, 1),
                      (0, 2),
                      (0, typemax(Int))]
        @test critical_edges == [DiameterSimplex{2}(3, 3, 1)]
    end

    @testset "sparse Float64" begin
        dist = sparse([0 2 0 0 5 0;
                       2 0 4 6 0 0;
                       0 4 0 3 0 0;
                       0 6 3 0 1 0;
                       5 0 0 1 0 0;
                       0 0 0 0 0 0])
        edges = copy_edges!(DiameterSimplex{3, Float64}[], dist, binomial)
        critical_edges = DiameterSimplex{3, Float64}[]

        res = compute_0_dim_pairs!(critical_edges, dist, binomial)

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
