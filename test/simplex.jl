using Ripserer: PrimeField, is_prime, n_bits, small_binomial

struct FakeFiltration <: AbstractFiltration{Int, Simplex{0,2,Int,UInt}} end
Ripserer.diam(::FakeFiltration, args...) =
    1
Ripserer.n_vertices(::FakeFiltration) =
    20

struct FakeFiltrationWithThreshold <: AbstractFiltration{Int, Simplex{0,2,Int,UInt}} end
Ripserer.diam(::FakeFiltrationWithThreshold, _, _, v) =
    v ≤ 10 ? 1 : ∞
Ripserer.n_vertices(::FakeFiltrationWithThreshold) =
    20

@testset "simplex.jl" begin
    @testset "helpers" begin
        @testset "is_prime" begin
            @test !is_prime(1)
            @test is_prime(2)
            @test is_prime(3)
            @test !is_prime(4)
            @test is_prime(5)
            @test !is_prime(6)
            @test is_prime(7)
            @test !is_prime(8)
            @test !is_prime(9)
            @test !is_prime(10)
            @test is_prime(11)
            @test !is_prime((1,2,3))
            @test !is_prime(:two)
            @test !is_prime(Array{Float64, 2})
        end
        @testset "n_bits" begin
            @test n_bits(2) == 1
            @test n_bits(3) == 2
            @test n_bits(7) == 3
            @test n_bits(Int64(typemax(UInt32))) == 32
        end
    end
    @testset "PrimeField" begin
        @test PrimeField{2}(1) + PrimeField{2}(1) == PrimeField{2}(0)
        @test PrimeField{5}(2) - PrimeField{5}(3) == PrimeField{5}(4)
        @test PrimeField{7}(2) / PrimeField{7}(3) == PrimeField{7}(3)
        @test PrimeField{17}(2) / PrimeField{17}(2) == PrimeField{17}(1)
        @test PrimeField{3}(2) * PrimeField{3}(2) == PrimeField{3}(1)

        @test PrimeField{3}(1) + 1 == PrimeField{3}(2)
        @test PrimeField{5}(2) * 2 == PrimeField{5}(4)
        @test PrimeField{7}(2) - 3 == PrimeField{7}(6)
        @test PrimeField{13}(10) / 2 == PrimeField{13}(5)

        for i in 1:10
            @test inv(PrimeField{11}(i)) == invmod(i, 11)
        end

        @test_throws DomainError PrimeField{4}(1)
        @test_throws DomainError PrimeField{-1}(1)

        @test sprint(print, PrimeField{3}(1)) == "1 mod 3"
    end
    @testset "Binomials" begin
        @test all(binomial(n, k) == small_binomial(n, Val(k)) for n in 0:1000 for k in 0:7)
    end
    @testset "Simplex" begin
        @testset "index, diam, coef, set_coef" begin
            for M in (2, 17, 7487), i in (1, 536, Int64(typemax(Int32)), Int32(10))
                c = rand(Int64)
                d = rand(Float64)
                @test index(Simplex{2, M}(d, i, c)) === i
                @test coef(Simplex{3, M}(d, i, c)) == PrimeField{M}(c)
                @test diam(Simplex{4, M}(d, i, c)) == d
                @test isa(Simplex{1, M}(d, i, c),
                          AbstractSimplex{1, PrimeField{M}, Float64})
            end
            @test_throws DomainError Simplex{2, -3}(rand(), 1, 1)
            @test_throws DomainError Simplex{2, 7497}(rand(), 1, 1)
            @test_throws DomainError Simplex{-1, 2}(rand(), 1, 1)

            @test set_coef(Simplex{1, 3}(1, 2, 1), 2) == Simplex{1, 3}(1, 2, 2)
            @test set_coef(Simplex{3, 2}(2, 2, 1), 3) == Simplex{3, 2}(2, 2, 1)

            for d in 1:10
                @test coface_type(Simplex{d, 3}(rand(), rand(Int), rand(Int))) ===
                    Simplex{d+1, 3, Float64, UInt}
                @test dim(Simplex{d, 2}(rand(), rand(Int), rand(Int))) == d
                @test dim(Simplex{d, 5, Float64, UInt64}) == d
            end

            @test sprint((io, val) -> show(io, MIME"text/plain"(), val),
                         Simplex{2,2}(1, 2, 3)) ==
                             """
                             2-dim Simplex{2}(1, 2, 1):
                               (4, 2, 1)"""
            @test sprint(print, Simplex{2,2}(1, 2, 3)) == "Simplex{2, 2}(1, (4, 2, 1), 1)"
        end
        @testset "arithmetic" begin
            @test Simplex{1, 3}(1.0, 3, 2) * 2 == Simplex{1, 3}(1.0, 3, 1)
            @test 2 * Simplex{2, 3}(2.0, 3, 2) == Simplex{2, 3}(2.0, 3, 1)
            @test -Simplex{3, 5}(10, 1, 1) == Simplex{3, 5}(10, 1, 4)

            @test Simplex{4, 7}(11, 2, 3) + Simplex{4, 7}(11, 2, 1) ==
                Simplex{4, 7}(11, 2, 4)

            @test Simplex{5, 2}(12, 3, 1) + Simplex{5, 2}(12, 3, 1) ==
                Simplex{5, 2}(12, 3, 0)

            for i in 1:16
                @test Simplex{2, 17}(i*10, 10, i) / i == Simplex{2, 17}(i*10, 10, 1)
                @test -Simplex{3, 17}(i*10, 8, i) == Simplex{3, 17}(i*10, 8, 17 - i)

                @test coef(Simplex{4, 17}(1.0, 1, i) - Simplex{4, 17}(1.0, 1, 17-i)) ==
                    coef(Simplex{4, 17}(1.0, 1, i) + -Simplex{4, 17}(1.0, 1, 17-i))
            end
        end
        @testset "vertices, index" begin
            @test vertices(Simplex{2, 2}(rand(Int), 1, 1)) == (3, 2, 1)
            @test vertices(Simplex{3, 2}(rand(Int), 2, 1)) == (5, 3, 2, 1)
            @test vertices(Simplex{1, 2}(rand(Int), 3, 1)) == (3, 2)
            @test vertices(Simplex{4, 2}(rand(Int), 4, 1)) == (6, 5, 4, 2, 1)
            @test vertices(Simplex{2, 2}(rand(Int), 5, 1)) == (5, 2, 1)

            for i in 1:20
                vxs = vertices(Simplex{5, 2}(rand(Int), i, 1))
                @test index(vxs) == i
            end
        end
        @testset "coboundary" begin
            @testset "number of cofaces" begin
                for dim in 1:10
                    cob = Simplex{dim+1, 2, Int, UInt}[]
                    simplex = Simplex{dim, 2}(1, 10, 1)
                    simplex_vxs = vertices(simplex)
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
                    cob = Simplex{dim+1, 2, Int, UInt}[]
                    for idx in binomial(20, dim+1):-1:1
                        simplex = Simplex{dim, 2}(1, idx, 1)
                        for coface in coboundary(FakeFiltration(), simplex, Val(false))
                            push!(cob, coface)
                        end
                    end
                    @test sort(index.(cob)) == 1:binomial(20, dim+2)
                end
            end
            @testset "number of cofaces, thresholding" begin
                for dim in 1:8 # at 9, simplex with index 10 becomes invalid.
                    cob = Simplex{dim+1, 2, Int}[]
                    simplex = Simplex{dim, 2}(1, 10, 1)
                    simplex_vxs = vertices(simplex)
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
                    simplex = Simplex{dim, 2}(1, 10, 1)
                    @inferred coboundary(FakeFiltration(), simplex)
                    cobiter = coboundary(FakeFiltration(), simplex)
                    @static if VERSION ≥ v"1.1.0"
                        @inferred Union{
                            Nothing,
                            Tuple{Simplex{dim+1, 2, Int, UInt}, Tuple{Int, Int}},
                        } iterate(cobiter)
                        @inferred Union{
                            Nothing,
                            Tuple{Simplex{dim+1, 2, Int, UInt}, Tuple{Int, Int}},
                        } iterate(cobiter, (13, dim+1))
                    end
                end
            end
        end
    end
end
