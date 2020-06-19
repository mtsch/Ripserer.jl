using Ripserer
using Ripserer: distances, vertex_type, edge_type, dist, edges, n_vertices

using Compat
using StaticArrays

using ..TestHelpers: test_filtration_interface
include("data.jl")

@testset "Edges for matrix inputs..." begin
    @testset "in a dense matrix." begin
        dist = [0 1 2 3;
                1 0 4 5;
                2 4 0 4;
                3 5 4 0]
        res = edges(dist, 4, Simplex{1, Int, Int})

        @test eltype(res) === Simplex{1, Int, Int}
        @test length(res) == 5
        @test res[1] == Simplex{1}(1, 1)
        @test res[2] == Simplex{1}(2, 2)
        @test res[3] == Simplex{1}(4, 3)
        @test res[4] == Simplex{1}(6, 4)
        @test res[5] == Simplex{1}(3, 4)
    end
    @testset "in a sparse matrix." begin
        dist = sparse([0 1 0 3;
                       1 0 4 5;
                       0 4 0 4;
                       3 5 4 0])
        res = edges(dist, Inf, Simplex{1, Int, Int128})

        @test eltype(res) === Simplex{1, Int, Int128}
        @test length(res) == 5
        @test res[1] === Simplex{1, Int, Int128}(1, 1)
        @test res[2] === Simplex{1, Int, Int128}(4, 3)
        @test res[3] === Simplex{1, Int, Int128}(6, 4)
        @test res[4] === Simplex{1, Int, Int128}(3, 4)
        @test res[5] === Simplex{1, Int, Int128}(5, 5)
    end
    @testset "number of edges in a dense matrix." begin
        dist = rand_dist_matrix(100)
        n_edges = binomial(size(dist, 1), 2)
        res = edges(dist, Inf, Simplex{1, Float64, Int})

        @test eltype(res) ≡ Simplex{1, Float64, Int}
        @test length(res) == n_edges
        @test issorted(res, by=diam)
    end
    @testset "number of edges in a sparse matrix." begin
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
    for points in (
        [(0, 0), (0, 1), (1, 1), (1, 0)],
        [SVector(0, 0), SVector(0, 1), SVector(1, 1), SVector(1, 0)]
    )
        @test distances(Euclidean(), points) ≈ [0 1 √2 1; 1 0 1 √2; √2 1 0 1; 1 √2 1 0]
        @test distances(Cityblock(), points) == [0 1 2 1; 1 0 1 2; 2 1 0 1; 1 2 1 0]
    end
end

for Filtration in (Rips, SparseRips)
    @testset "$(string(Filtration))" begin
        test_filtration_interface(Filtration, (icosahedron, torus(100)))

        @testset "With no threshold." begin
            filtration = Filtration(Float64[0 1 2 9; 1 0 3 9; 2 3 0 4; 9 9 4 0])

            @test n_vertices(filtration) == 4
            @test dist(filtration, 3, 3) == 0.0
            @test dist(filtration, 1, 2) == 1.0
            @test dist(filtration, 1, 3) == 2.0
            @test dist(filtration, 3, 2) == 3.0
            @test diam(filtration, Simplex{1}(3, 3), [3, 2], 1) == 3.0
            @test diam(filtration, (4, 2, 1)) === (filtration isa Rips ? missing : 9.0)

            @test threshold(filtration) == (filtration isa Rips ? 4.0 : 9.0)
            @test vertex_type(filtration) === Simplex{0, Float64, Int}
            @test edge_type(filtration) === Simplex{1, Float64, Int}
            @test birth(filtration, 1) == 0
            @test birth(filtration, 4) == 0
        end
        @testset "With threshold and index type." begin
            filtration = Filtration{Int128}([1 1 2; 1 2 3; 2 3 2]; threshold=2)

            @test n_vertices(filtration) == 3
            @test threshold(filtration) == 2
            @test dist(filtration, 3, 3) == 0
            @test dist(filtration, 1, 2) == 1
            @test dist(filtration, 1, 3) == 2
            @test dist(filtration, 3, 2) ≡ (issparse(filtration.dist) ? missing : 3)
            @test ismissing(diam(filtration, Simplex{2}(1, 1), [1, 2], 3))
            @test threshold(filtration) == 2
            @test vertex_type(filtration) === Simplex{0, Int, Int128}
            @test edge_type(filtration) === Simplex{1, Int, Int128}

            @test birth(filtration, 1) == 1
            @test birth(filtration, 2) == 2
            @test birth(filtration, 3) == 2
        end
    end
end

@testset "Rips points constructor." begin
    filtration = Rips(torus_points(9))
    @test all(x -> x > 0, dist(filtration, i, j) for i in 1:10 for j in i+1:9)
    @test eltype(edges(filtration)) === Simplex{1, Float64, Int}

    filtration = Rips{Int32}(torus_points(9))
    @test all(x -> x > 0, dist(filtration, i, j) for i in 1:10 for j in i+1:9)
    @test eltype(edges(filtration)) === Simplex{1, Float64, Int32}
end

@testset "Errors" begin
    @testset "Non-distance matrices throw an error." begin
        @test_throws ArgumentError Rips([1 2 3; 4 5 6; 7 8 9])
        @test_throws ArgumentError Rips(zeros(3, 2))
        @test_throws ArgumentError SparseRips([1 2 3; 4 5 6; 7 8 9])
        @test_throws ArgumentError SparseRips(zeros(3, 2))
    end
    @testset "Constructing Rips filtration with sparse matrix not allowed." begin
        @test_throws ArgumentError Rips(sparse([0 1 1; 1 0 1; 1 1 0]))
    end
end
