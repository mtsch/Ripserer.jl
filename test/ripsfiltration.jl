using Ripserer
using Ripserer: distances, dist_type, vertex_type, edge_type, dist, threshold

include("data.jl")

@testset "edges" begin
    @testset "edges dense" begin
        dist = [0 1 2 3;
                1 0 4 5;
                2 4 0 4;
                3 5 4 0]
        res = edges(dist, ∞, Simplex{1, 2, Int, UInt64})
        @test eltype(res) isa DataType
        @test res[1] == Simplex{1, 2}(1, 1, 1)
        @test res[2] == Simplex{1, 2}(2, 2, 1)
        @test res[3] == Simplex{1, 2}(3, 4, 1)
        @test res[4] == Simplex{1, 2}(4, 6, 1)
        @test res[5] == Simplex{1, 2}(4, 3, 1)
        @test res[6] == Simplex{1, 2}(5, 5, 1)

        res = edges(dist, 3, Simplex{1, 2, Int, UInt64})
        @test eltype(res) isa DataType
        @test res[1] == Simplex{1, 2}(1, 1, 1)
        @test res[2] == Simplex{1, 2}(2, 2, 1)
        @test res[3] == Simplex{1, 2}(3, 4, 1)
    end
    @testset "edges sparse" begin
        dist = sparse([0 1 0 3;
                       1 0 4 5;
                       0 4 0 4;
                       3 5 4 0])
        res = edges(dist, ∞, Simplex{1, 3, Int, UInt64})
        @test eltype(res) isa DataType
        @test res[1] == Simplex{1, 3}(1, 1, 1)
        @test res[2] == Simplex{1, 3}(3, 4, 1)
        @test res[3] == Simplex{1, 3}(4, 6, 1)
        @test res[4] == Simplex{1, 3}(4, 3, 1)
        @test res[5] == Simplex{1, 3}(5, 5, 1)
    end
    @testset "n edges dense" begin
        dist = rand_dist_matrix(100)
        n_edges = binomial(size(dist, 1), 2)
        @test length(edges(dist, ∞, Simplex{1, 7})) == n_edges
    end
    @testset "n edges sparse" begin
        for _ in 1:10
            dist = rand_dist_matrix(100, 0.5)
            n_edges = nnz(dist) ÷ 2
            @test length(edges(dist, ∞, Simplex{1, 2, Float64, UInt})) == n_edges
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

"""
    filtration_interface_testset(Filtration, datasets)

Test the filtration interface. Make sure all required methods are defined and that they
return values of types that make sense.
"""
function test_filtration_interface(Filtration, datasets)
    @testset "AbstractFiltration interface" begin
        for data in datasets
            flt = Filtration(data)

            T = dist_type(flt)
            V = vertex_type(flt)
            E = edge_type(flt)

            @test V <: AbstractSimplex{0}
            @test V isa DataType
            @test E <: AbstractSimplex{1}
            @test E isa DataType

            @test n_vertices(flt) == size(data, 1)
            sx = first(edges(flt))
            @test sx isa E
            @test diam(sx) isa T

            @test diam(flt, [3, 2, 1]) isa T
            @test diam(flt, E(one(T), 1, 1), [3, 2], 1) isa T
            @test issparse(flt) isa Bool
            @test issparse(Filtration) isa Bool
            @test issparse(Filtration) == issparse(flt)
        end
    end
end

for Filtration in (RipsFiltration, SparseRipsFiltration)
    @testset "$(string(Filtration))" begin
        test_filtration_interface(Filtration, (icosahedron, torus(100)))

        @testset "no threshold, modulus=3" begin
            flt = Filtration([0 1 2 9;
                              1 0 3 9;
                              2 3 0 4;
                              9 9 4 0], modulus=3)
            @test n_vertices(flt) == 4
            @test dist(flt, 3, 3) == 0
            @test dist(flt, 1, 2) == 1
            @test dist(flt, 1, 3) == 2
            @test dist(flt, 3, 2) == 3
            @test diam(flt, Simplex{1, 3}(3, 3, 1), [3, 2], 1) == 3
            @test threshold(flt) == (issparse(Filtration) ? ∞ : 4)
            @test vertex_type(flt) === Simplex{0, 3, Int, UInt}
            @test edge_type(flt) === Simplex{1, 3, Int, UInt}
            @test birth(flt, 1) == 0
            @test birth(flt, 4) == 0
        end
        @testset "threshold, vertex_type" begin
            flt = Filtration([1 1 2;
                              1 2 3;
                              2 3 3];
                             threshold=2,
                             vertex_type=Simplex{0, 5, Int, UInt})

            @test n_vertices(flt) == 3
            @test dist(flt, 3, 3) == 0
            @test dist(flt, 1, 2) == 1
            @test dist(flt, 1, 3) == 2
            @test dist(flt, 3, 2) == (issparse(Filtration) ? ∞ : 3)
            @test diam(flt, Simplex{2, 5}(1, 1, 1), [1, 2], 3) == ∞
            @test threshold(flt) == (issparse(Filtration) ? ∞ : 2)
            @test vertex_type(flt) === Simplex{0, 5, Int, UInt}
            @test edge_type(flt) === Simplex{1, 5, Int, UInt}

            @test birth(flt, 1) == 1
            @test birth(flt, 2) == 2
            @test birth(flt, 3) == 3
        end
        @testset "errors" begin
            @test_throws ArgumentError Filtration([1 2 3;
                                                   4 5 6;
                                                   7 8 9])
            dists = [0 1 2; 1 0 3; 2 3 0]
            @test_throws ArgumentError Filtration(dists, modulus=4)
            @test_throws ArgumentError Filtration(dists, modulus=-2)
            @test_throws ArgumentError Filtration(dists, vertex_type=Int)
            @test_throws ArgumentError Filtration(dists, vertex_type=Simplex{1,2,Int,UInt})
            @test_throws TypeError Filtration(dists, vertex_type=Simplex{0, 2, Int})
        end
    end
end
