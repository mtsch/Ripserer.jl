using Ripserer
using Ripserer: small_binomial

using ..TestHelpers: test_indexed_simplex_interface

struct FakeFiltration end
Ripserer.diam(::FakeFiltration, args...) =
    1
Ripserer.n_vertices(::FakeFiltration) =
    20

struct FakeFiltrationWithThreshold end
Ripserer.diam(::FakeFiltrationWithThreshold, _, _, v) =
    v ≤ 10 ? 1 : ∞
Ripserer.n_vertices(::FakeFiltrationWithThreshold) =
    20

@testset "Binomials" begin
    @test all(binomial(n, k) == small_binomial(n, Val(k)) for n in 0:1000 for k in 0:7)
end

@testset "Simplex" begin
    test_indexed_simplex_interface(Simplex, D->D+1)

    @testset "vertices, index" begin
        @test vertices(Simplex{2}(1, rand())) == (3, 2, 1)
        @test vertices(Simplex{3}(2, rand())) == (5, 3, 2, 1)
        @test vertices(Simplex{1}(3, rand())) == (3, 2)
        @test vertices(Simplex{4}(4, rand())) == (6, 5, 4, 2, 1)
        @test vertices(Simplex{2}(5, rand())) == (5, 2, 1)

        for i in 1:20
            sx = Simplex{5}(i, rand())
            @test index(vertices(sx)) == i
            @test Simplex{5}(vertices(sx), diam(sx)) == sx
            @test_throws ArgumentError Simplex{4}(vertices(sx), rand())
            @test_throws ArgumentError Simplex{6}(vertices(sx), rand())
        end
    end
    @testset "show" begin
        @test sprint(print, Simplex{1}(1, 1)) == "Simplex{1}(+(2, 1), 1)"
        @test sprint(print, Simplex{2}(-1, 1)) == "Simplex{2}(-(3, 2, 1), 1)"

        @test sprint(Simplex{2}(1, 1)) do io, sx
            show(io, MIME"text/plain"(), sx)
        end == "2-dim Simplex(1, 1):\n  +(3, 2, 1)"

        @test sprint(Simplex{1}(Int128(1), 1)) do io, sx
            show(io, MIME"text/plain"(), sx)
        end == "1-dim Simplex(1, 1) with Int128 index:\n  +(2, 1)"

    end
    @testset "coboundary" begin
        @testset "number of cofaces" begin
            for dim in 1:10
                simplex = Simplex{dim}(10, 1)
                simplex_vxs = vertices(simplex)
                cob = coface_type(simplex)[]
                for coface in coboundary(FakeFiltration(), simplex)
                    push!(cob, coface)
                end
                @test length(cob) == 20 - dim - 1
                @test all(isequal(1), diam.(cob))
                @test all(issubset(simplex_vxs, vertices(c)) for c in cob)
            end
        end
        @testset "number of cofaces, all_cofaces=false" begin
            for dim in 1:5
                cob = Simplex{dim+1, Int, Int}[]
                for idx in binomial(20, dim+1):-1:1
                    simplex = Simplex{dim}(idx, 1)
                    for coface in coboundary(FakeFiltration(), simplex, Val(false))
                        push!(cob, coface)
                    end
                end
                @test sort(index.(abs.(cob))) == 1:binomial(20, dim+2)
            end
        end
        @testset "number of cofaces, thresholding" begin
            for dim in 1:8 # at 9, simplex with index 10 becomes invalid.
                simplex = Simplex{dim}(10, 1)
                simplex_vxs = vertices(simplex)
                cob = coface_type(simplex)[]
                for coface in coboundary(FakeFiltrationWithThreshold(), simplex)
                    push!(cob, coface)
                end
                @test length(cob) == 10 - dim - 1
                @test all(isequal(1), diam.(cob))
                @test all(issubset(simplex_vxs, vertices(c)) for c in cob)
            end
        end
        @testset "type stability" begin
            for dim in 1:10
                simplex = Simplex{dim}(10, 1)
                @test begin @inferred coboundary(FakeFiltration(), simplex); true end
                cobiter = coboundary(FakeFiltration(), simplex)
                @static if VERSION ≥ v"1.1.0"
                    @test begin @inferred Union{
                        Nothing,
                        Tuple{coface_type(simplex), Tuple{Int, Int}},
                    } iterate(cobiter); true end
                    @test begin @inferred Union{
                        Nothing,
                        Tuple{coface_type(simplex), Tuple{Int, Int}},
                    } iterate(cobiter, (13, dim+1)); true end
                end
            end
        end
    end
end