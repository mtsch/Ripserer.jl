using Ripserer
using Ripserer: small_binomial

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
    @testset "interface" begin
        for D in 0:5
            for T in (Float64, Float32, Int)
                d = rand(T)
                for I in (Int64, Int128)
                    i = rand(I)

                    @test index(Simplex{D}(i, d)) == i
                    @test diam(Simplex{D}(i, d)) == d
                    @test typeof(Simplex{D}(i, d)) ≡ Simplex{D, T, I}
                    @test Simplex{D}(i, d) isa IndexedSimplex{D, T, I}
                    @test Simplex{D}(i, d) isa AbstractSimplex{D, T}
                    D > 0 && @test_throws DomainError Simplex{-D}(i, d)

                    @test -Simplex{D}(i, d) == Simplex{D}(-i, d)
                    @test sign(+Simplex{D}(i, d)) == sign(i)
                    @test sign(-Simplex{D}(i, d)) == -sign(i)

                    @test coface_type(Simplex{D}(i, d)) ≡ Simplex{D + 1, T, I}
                    @test coface_type(typeof(Simplex{D}(i, d))) ≡ Simplex{D + 1, T, I}

                    @test dim(Simplex{D}(i, d)) == D
                    @test abs(Simplex{D}(i, d)) == Simplex{D}(abs(i), d)
                    @test abs(-Simplex{D}(i, d)) == abs(Simplex{D}(i, d))
                end
                # don't want to do this with random int value.
                @test length(vertices(Simplex{D}(10, d))) == D + 1
                @test length(vertices(Simplex{D}(-10, d))) == D + 1

                @inferred vertices(Simplex{D}(10, d))
                @inferred vertices(Simplex{D}(-10, d))
            end
        end

        for sgn in (1, -1)
            @test isless(Simplex{0}(sgn * 10, 1), Simplex{0}(10, 2))
            @test !isless(Simplex{1}(sgn * 10, 2), Simplex{1}(10, 1))
            @test isless(Simplex{2}(sgn * 11, 1), Simplex{2}(10, 1))
            @test !isless(Simplex{3}(sgn * 9, 1), Simplex{3}(10, 1))
            @test !isless(Simplex{4}(sgn * 10, 1), Simplex{4}(10, 1))
        end
    end
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
