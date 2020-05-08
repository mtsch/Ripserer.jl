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

@testset "Mod" begin
    @test Mod{2}(1) + Mod{2}(1) == Mod{2}(0)
    @test Mod{5}(2) - Mod{5}(3) == Mod{5}(4)
    @test Mod{7}(2) / Mod{7}(3) == Mod{7}(3)
    @test Mod{17}(2) / Mod{17}(2) == Mod{17}(1)
    @test Mod{3}(2) * Mod{3}(2) == Mod{3}(1)

    @test Mod{3}(1) + 1 == Mod{3}(2)
    @test Mod{5}(2) * 2 == Mod{5}(4)
    @test Mod{7}(2) - 3 == Mod{7}(6)
    @test Mod{13}(10) / 2 == Mod{13}(5)

    for i in 1:10
        @test inv(Mod{11}(i)) == invmod(i, 11)
    end

    @test_throws DomainError Mod{4}(1)
    @test_throws DomainError Mod{-1}(1)

    @test sprint(print, Mod{3}(1)) == "1 mod 3"
end
