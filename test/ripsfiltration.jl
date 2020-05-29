using Ripserer
using Ripserer: distances, dist_type, vertex_type, edge_type, dist, threshold

using ..TestHelpers: test_filtration_interface
include("data.jl")

@testset "edges" begin
    @testset "edges dense" begin
        dist = [0 1 2 3;
                1 0 4 5;
                2 4 0 4;
                3 5 4 0]
        res = edges(dist, Inf, Simplex{1, Int, Int})
        @test eltype(res) ≡ Simplex{1, Int, Int}
        @test res[1] == Simplex{1}(1, 1)
        @test res[2] == Simplex{1}(2, 2)
        @test res[3] == Simplex{1}(4, 3)
        @test res[4] == Simplex{1}(6, 4)
        @test res[5] == Simplex{1}(3, 4)
        @test res[6] == Simplex{1}(5, 5)

        res = edges(dist, 3, Simplex{1})
        @test eltype(res) isa DataType
        @test res[1] == Simplex{1}(1, 1)
        @test res[2] == Simplex{1}(2, 2)
        @test res[3] == Simplex{1}(4, 3)
    end
    @testset "edges sparse" begin
        dist = sparse([0 1 0 3;
                       1 0 4 5;
                       0 4 0 4;
                       3 5 4 0])
        res = edges(dist, Inf, Simplex{1, Int, Int128})
        @test eltype(res) ≡ Simplex{1, Int, Int128}
        @test res[1] == Simplex{1}(1, 1)
        @test res[2] == Simplex{1}(4, 3)
        @test res[3] == Simplex{1}(6, 4)
        @test res[4] == Simplex{1}(3, 4)
        @test res[5] == Simplex{1}(5, 5)
    end
    @testset "n edges dense" begin
        dist = rand_dist_matrix(100)
        n_edges = binomial(size(dist, 1), 2)
        res = edges(dist, Inf, Simplex{1})
        @test eltype(res) ≡ Simplex{1, Float64, Int}
        @test length(res) == n_edges
        @test issorted(res, by=diam)
    end
    @testset "n edges sparse" begin
        for _ in 1:10
            dist = rand_dist_matrix(100, 0.5)
            n_edges = nnz(dist) ÷ 2
            res = edges(dist, Inf, Simplex{1, Float64, Int128})
            @test eltype(res) ≡ Simplex{1, Float64, Int128}
            @test length(res) == n_edges
            @test issorted(res, by=diam)
        end
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

for Filtration in (Rips, SparseRips)
    @testset "$(string(Filtration))" begin
        test_filtration_interface(Filtration, (icosahedron, torus(100)))

        @testset "no threshold" begin
            flt = Filtration(Float64[0 1 2 9;
                                     1 0 3 9;
                                     2 3 0 4;
                                     9 9 4 0])
            @test n_vertices(flt) == 4
            @test dist(flt, 3, 3) == 0.0
            @test dist(flt, 1, 2) == 1.0
            @test dist(flt, 1, 3) == 2.0
            @test dist(flt, 3, 2) == 3.0
            @test diam(flt, Simplex{1}(3, 3), [3, 2], 1) == 3.0
            @test threshold(flt) == (issparse(flt.dist) ? 9 : 4.0)
            @test vertex_type(flt) === Simplex{0, Float64, Int}
            @test edge_type(flt) === Simplex{1, Float64, Int}
            @test birth(flt, 1) == 0
            @test birth(flt, 4) == 0
        end
        @testset "threshold, vertex_type" begin
            flt = Filtration([1 1 2;
                              1 2 3;
                              2 3 3];
                             threshold=2,
                             vertex_type=Simplex{0, Int, Int})

            @test n_vertices(flt) == 3
            @test threshold(flt) == 2
            @test dist(flt, 3, 3) == 0
            @test dist(flt, 1, 2) == 1
            @test dist(flt, 1, 3) == 2
            @test dist(flt, 3, 2) ≡ (issparse(flt.dist) ? missing : 3)
            @test ismissing(diam(flt, Simplex{2}(1, 1), [1, 2], 3))
            @test threshold(flt) == 2
            @test vertex_type(flt) === Simplex{0, Int, Int}
            @test edge_type(flt) === Simplex{1, Int, Int}

            @test birth(flt, 1) == 1
            @test birth(flt, 2) == 2
            @test birth(flt, 3) == 3
        end
        @testset "errors" begin
            @test_throws ArgumentError Filtration([1 2 3;
                                                   4 5 6;
                                                   7 8 9])
            dists = [0 1 2; 1 0 3; 2 3 0]
            @test_throws ArgumentError Filtration(dists, vertex_type=Int)
            @test_throws ArgumentError Filtration(
                dists, vertex_type=Simplex{1, Int, Int}
            )
            @test_throws ArgumentError Filtration(
                dists, vertex_type=Simplex{0, Float64, Int}
            )
            @test_throws TypeError Filtration(dists, vertex_type=Simplex{0, Int})
        end
    end
end
