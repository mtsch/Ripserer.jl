using Ripserer
using Ripserer: is_prime

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
