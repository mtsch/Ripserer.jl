using Ripserer: PrimeField, isprime, n_bits, set_coef

@testset "simplex.jl" begin
    @testset "helpers" begin
        @testset "isprime" begin
            @test !isprime(1)
            @test isprime(2)
            @test isprime(3)
            @test !isprime(4)
            @test isprime(5)
            @test !isprime(6)
            @test isprime(7)
            @test !isprime(8)
            @test !isprime(9)
            @test !isprime(10)
            @test isprime(11)
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
    end

    @testset "Simplex" begin
        @testset "index, diam, coef, set_coef" begin
            for M in (2, 17, 7487)
                for i in (1, 536, typemax(Int32))
                    c = rand(Int64)
                    d = rand(Float64)
                    @test Simplex{M}(d, i, c) isa AbstractSimplex{PrimeField{M}, Float64}
                    @test index(Simplex{M}(d, i, c)) == i
                    @test coef(Simplex{M}(d, i, c)) == PrimeField{M}(c)
                    @test diam(Simplex{M}(d, i, c)) == d
                end
            end
            @test_throws DomainError Simplex{-3}(rand(), 1, 1)
            @test_throws DomainError Simplex{7497}(rand(), 1, 1)

            @test set_coef(Simplex{3}(1, 2, 1), 2) == Simplex{3}(1, 2, 2)
            @test set_coef(Simplex{2}(2, 2, 1), 3) == Simplex{2}(2, 2, 1)
        end

        @testset "arithmetic" begin
            @test Simplex{3}(1.0, 3, 2) * 2 == Simplex{3}(1.0, 3, 1)
            @test 2 * Simplex{3}(2.0, 3, 2) == Simplex{3}(2.0, 3, 1)
            @test -Simplex{5}(10, 1, 1) == Simplex{5}(10, 1, 4)
            @test Simplex{7}(11, 2, 3) + Simplex{7}(11, 2, 1) == Simplex{7}(11, 2, 4)
            @test Simplex{2}(12, 3, 1) + Simplex{2}(12, 3, 1) == Simplex{2}(12, 3, 0)

            for i in 1:16
                @test Simplex{17}(i*10, 10, i) / i == Simplex{17}(i*10, 10, 1)
                @test -Simplex{17}(i*10, 8, i) == Simplex{17}(i*10, 8, 17 - i)
                @test coef(Simplex{17}(1.0, 1, i) - Simplex{17}(1.0, 1, 17-i)) ==
                    coef(Simplex{17}(1.0, 1, i) + -Simplex{17}(1.0, 1, 17-i))
            end
        end
    end
end
