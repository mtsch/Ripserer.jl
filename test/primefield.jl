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
    @testset "Arithmetic" begin
        @test Mod{2}(1) + Mod{2}(1) == Mod{2}(0)
        @test Mod{5}(2) - Mod{5}(3) == Mod{5}(4)
        @test Mod{7}(2) / Mod{7}(3) == Mod{7}(3)
        @test Mod{17}(2) / Mod{17}(2) == Mod{17}(1)
        @test Mod{3}(2) * Mod{3}(2) == Mod{3}(1)

        @test Mod{3}(1) + 1 == Mod{3}(2)
        @test Mod{5}(2) * 2 == Mod{5}(4)
        @test Mod{7}(2) - 3 == Mod{7}(6)
        @test Mod{13}(10) / 2 == Mod{13}(5)

        @test sign(Mod{2}(1)) == Mod{2}(1)
        @test sign(Mod{13}(0)) == Mod{13}(0)

        @test zero(Mod{2}) == Mod{2}(0)
        @test zero(Mod{5}(1)) == Mod{5}(0)

        for i in 1:10
            @test inv(Mod{11}(i)) == invmod(i, 11)
        end
    end

    @testset "Errors" begin
        @test_throws DomainError Mod{4}(1)
        @test_throws DomainError Mod{-1}(1)
        @test_throws DomainError Mod{String}(1)
        @test_throws DomainError Mod{-1}(1)

        @test_throws ErrorException Mod{3}(1) < Mod{3}(2)
        @test_throws ErrorException Mod{3}(1) ≥ Mod{3}(2)

        @test_throws StackOverflowError Mod{3}(1) + Mod{2}(1)
    end

    @testset "Printing" begin
        @test sprint(print, Mod{3}(1)) == "1 mod 3"
        @test sprint(print, Mod{1801}(1)) == "1 mod 1801"
    end

    @testset "Type promotion" begin
        for T in (Int8, Int16, Int32, Int64, Int128)
            @test promote_type(T, Mod{2}) ≡ Mod{2}
            @test promote(one(T), Mod{5}(2)) == (Mod{5}(1), Mod{5}(2))

            @test Mod{3}(one(T)) == Mod{3}(1)
            @test Mod{3}(one(unsigned(T))) == Mod{3}(1)
            for op in (+, -, *, /)
                @test op(Mod{3}(1), one(T)) isa Mod{3}
                @test op(Mod{3}(1), one(unsigned(T))) isa Mod{3}
            end
        end
    end
end
