using Ripserer: isprime, Binomial

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

    @testset "Binomial" begin
        bin = Binomial(10, 15)
        for n in 1:10, k in 1:15
            @test bin(n, k) == binomial(n, k)
        end
    end
end
