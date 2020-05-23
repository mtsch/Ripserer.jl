using Ripserer
using Ripserer: adj_matrix,
    BottleneckGraph, depth_layers, augmenting_paths, augment!, hopcroft_karp!

@testset "adj_matrix" begin
    diag1 = PersistenceDiagram(0, [(1, 2), (1, ∞)])
    diag2 = PersistenceDiagram(0, [(3, 4), (5, 10), (7, ∞)])

    #                                  1,2 1,∞ 3,3 5,5 7,7
    @test adj_matrix(diag1, diag2) == [2.0 Inf 1.0 Inf Inf  # 3,4
                                       8.0 Inf Inf 5.0 Inf  # 5,7
                                       Inf 6.0 Inf Inf Inf  # 7,∞
                                       1.0 Inf 0.0 0.0 0.0  # 1,1
                                       Inf Inf 0.0 0.0 0.0] # 1,1

    @test adj_matrix(diag1, diag2) == adj_matrix(diag2, diag1)'
end

@testset "Hopcroft-Karp" begin
    adj = [1 1 9;
           1 9 1;
           9 9 2]
    graph = BottleneckGraph(adj, [0, 0, 0], [0, 0, 0], Int[], 3)

    @test depth_layers(graph, 2) == ([1, 1, 1], 1)
    @test augmenting_paths(graph, 2) == [[1, 1], [3, 2]]

    graph = BottleneckGraph(adj, [1, 0, 2], [1, 3, 0], Int[], 3)
    @test depth_layers(graph, 2) == ([1, 2, 3], 3)
    @test augmenting_paths(graph, 2) == [[2, 1, 1, 2, 3, 3]]
    augment!(graph, [2, 1, 1, 2, 3, 3])
    @test graph.match_left == [2, 1, 3]
    @test graph.match_right == [2, 1, 3]

    graph = BottleneckGraph(adj, [0, 0, 0], [0, 0, 0], Int[], 3)
    @test hopcroft_karp!(graph, 2) == ([(1, 2), (2, 1), (3, 3)], true)
    @test hopcroft_karp!(graph, 1) == ([(1, 1), (3, 2)], false)
end

@testset "Bottleneck basic" begin
    diag1 = PersistenceDiagram(0, [(1, 2), (5, 8)])
    diag2 = PersistenceDiagram(0, [(1, 2), (3, 4), (5, 10)])

    m = matching(Bottleneck(), diag1, diag2)
    @test matching(m) == [(5,8) => (5,10)]
    @test matching(m, bottleneck=false) == [(1,2) => (1,2), (3,3) => (3,4), (5,8) => (5,10)]
    @test distance(m) ≡ 2
    @test distance(Bottleneck(), diag1, diag2) ≡ 2
    @test distance(Bottleneck(), diag1, diag2) == distance(Bottleneck(), diag2, diag1)

    @test distance(Bottleneck(), diag1, diag1) ≡ 0
end

@testset "Wasserstein basic" begin
    diag1 = PersistenceDiagram(0, [(1, 2), (5, 8)])
    diag2 = PersistenceDiagram(0, [(1, 2), (3, 4), (5, 10)])

    m = matching(Wasserstein(), diag1, diag2)
    @test matching(m) == [(1,2) => (1,2), (3,3) => (3,4), (5,8) => (5,10)]
    @test distance(m) ≡ 3
    @test distance(Wasserstein(), diag1, diag2) ≡ 3
    @test distance(Wasserstein(2), diag1, diag2) ≡ √(1 + 4)
    for i in 1:3
        @test distance(Wasserstein(i), diag1, diag2) ≡
            distance(Wasserstein(i), diag2, diag1)
    end

    @test distance(Wasserstein(), diag1, diag1) ≡ 0
end

@testset "infinite intervals" begin
    diag1 = PersistenceDiagram(0, [(1, 2), (5, 8), (1, ∞)])
    diag2 = PersistenceDiagram(0, [(1, 2), (3, 4), (5, 10)])
    diag3 = PersistenceDiagram(0, [(1, 2), (3, 4), (5, 10), (1, ∞)])
    diag4 = PersistenceDiagram(1, [(1, ∞)])
    diag5 = PersistenceDiagram(1, [(2, ∞)])

    for dist_type in (Bottleneck(), Wasserstein(), Wasserstein(2))
        @test distance(dist_type, diag1, diag2) ≡ ∞
        @test distance(dist_type, diag2, diag1) ≡ ∞
        @test distance(dist_type, diag1, diag1) == 0
        @test 0 < distance(dist_type, diag1, diag3) < ∞
        @test distance(dist_type, diag4, diag5) == 1
    end

end

@testset "different sizes" begin
    diag1 = PersistenceDiagram(5, vcat((90, 100), [(i, i+1) for i in 1:100]))
    diag2 = PersistenceDiagram(5, [(100, 110)])

    @test distance(Bottleneck(), diag1, diag2) == 10
    @test distance(Bottleneck(), diag2, diag1) == 10
    @test distance(Wasserstein(), diag1, diag2) == 110
    @test distance(Wasserstein(), diag2, diag1) == 110
    @test distance(Wasserstein(2), diag1, diag2) == √200
    @test distance(Wasserstein(2), diag2, diag1) == √200
end
