using Ripserer

@testset "Bottleneck" begin
    diag1 = PersistenceDiagram(0, [(1, 2), (5, 8)])
    diag2 = PersistenceDiagram(0, [(1, 2), (3, 4), (5, 10)])

    m = matching(Bottleneck(), diag1, diag2)
    @test matching(m) == [(1,2) => (1,2), (3,3) => (3,4), (5,8) => (5,10)]
    @test distance(m) ≡ 2
    @test distance(Bottleneck(), diag1, diag2) ≡ 2

    @test distance(Bottleneck(), diag1, diag1) ≡ 0
end

@testset "Wasserstein" begin
    diag1 = PersistenceDiagram(0, [(1, 2), (5, 8)])
    diag2 = PersistenceDiagram(0, [(1, 2), (3, 4), (5, 10)])

    m = matching(Wasserstein(), diag1, diag2)
    @test matching(m) == [(1,2) => (1,2), (3,3) => (3,4), (5,8) => (5,10)]
    @test distance(m) ≡ 3
    @test distance(Wasserstein(), diag1, diag2) ≡ 3
    @test distance(Wasserstein(2), diag1, diag2) ≡ √(1 + 4)

    @test distance(Wasserstein(), diag1, diag1) ≡ 0
end
