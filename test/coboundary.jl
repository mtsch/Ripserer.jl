using Ripserer: Binomial, Coboundary, vertices

@testset "coboundary" begin
    @testset "Binomial" begin
        bin = Binomial(10, 15)
        for n in 0:10, k in 0:15
            @test bin(n, k) == binomial(n, k)
        end
        @test_throws BoundsError bin(11, 15)
        @test_throws BoundsError bin(10, 16)
        @test_throws BoundsError bin(-1, 15)
        @test_throws BoundsError bin(10, -1)
    end

    @testset "coboundary" begin
        struct FakeFiltration <: AbstractFiltration{Int, Simplex{2, Int}} end
        Ripserer.diam(::FakeFiltration, args...) =
            1
        Ripserer.n_vertices(::FakeFiltration) =
            20

        struct FakeFiltrationWithThreshold <: AbstractFiltration{Int, Simplex{2, Int}} end
        Ripserer.diam(::FakeFiltrationWithThreshold, _, _, v) =
            v ≤ 10 ? 1 : ∞
        Ripserer.n_vertices(::FakeFiltrationWithThreshold) =
            20

        @testset "vertices, index" begin
            coboundary = Coboundary(FakeFiltration(), 10)
            @test vertices(coboundary, Simplex{2}(rand(Int), 1, 1), 2) == [3, 2, 1]
            @test vertices(coboundary, Simplex{2}(rand(Int), 2, 1), 3) == [5, 3, 2, 1]
            @test vertices(coboundary, Simplex{2}(rand(Int), 3, 1), 1) == [3, 2]
            @test vertices(coboundary, Simplex{2}(rand(Int), 4, 1), 4) == [6, 5, 4, 2, 1]
            @test vertices(coboundary, Simplex{2}(rand(Int), 5, 1), 2) == [5, 2, 1]

            for i in 1:10
                vxs = vertices(coboundary, Simplex{2}(rand(Int), i, 1), 5)
                @test index(coboundary, vxs) == i
            end
        end
        @testset "number of cofaces" begin
            coboundary = Coboundary(FakeFiltration(), 10)

            for dim in 1:10
                cob = Simplex{2, Int}[]
                simplex = Simplex{2}(1, 10, 1)
                simplex_vxs = vertices(coboundary, simplex, dim)
                for coface in coboundary(Simplex{2}(1, 1, 1), dim)
                    push!(cob, coface)
                end
                @test length(cob) == 20 - dim - 1
                @test all(isequal(1), diam.(cob))
                @test all(issubset(simplex_vxs, vertices(coboundary, c, dim+1))
                          for c in cob)
            end
        end
        @testset "number of cofaces, all_cofaces=false" begin
            coboundary = Coboundary(FakeFiltration(), 5)

            for dim in 1:5
                cob = Simplex{2, Int}[]
                for idx in binomial(20, dim+1):-1:1
                    simplex = Simplex{2}(1, idx, 1)
                    for coface in coboundary(simplex, dim, Val(false))
                        push!(cob, coface)
                    end
                end
                @test length(cob) == binomial(20, dim+2)
                @test allunique(cob)
            end
        end
        @testset "number of cofaces, thresholding" begin
            coboundary = Coboundary(FakeFiltrationWithThreshold(), 10)

            for dim in 1:10
                cob = Simplex{2, Int}[]
                simplex = Simplex{2}(1, 10, 1)
                simplex_vxs = vertices(coboundary, simplex, dim)
                for coface in coboundary(Simplex{2}(1, 1, 1), dim)
                    push!(cob, coface)
                end
                @test length(cob) == max(0, 10 - dim - 1)
                @test all(isequal(1), diam.(cob))
                @test all(issubset(simplex_vxs, vertices(coboundary, c, dim+1))
                          for c in cob)
            end
        end
    end
end
