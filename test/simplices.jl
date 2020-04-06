using Ripserer: Simplex, DiameterSimplex, index, coef, set_coef, diam, vertices, inv_mod

@testset "simplices" begin
    @testset "Simplex" begin
        for M in (2, 17, 7487)
            for i in (1, 536, typemax(Int32))
                c = rand(Int64)
                @test index(Simplex{M}(i, c)) == i
                @test coef(Simplex{M}(i, c)) == mod(c, M)
            end
        end
        @test_throws DomainError Simplex{-3}(1, 1)
        @test_throws DomainError Simplex{7497}(1, 1)
    end

    @testset "DiameterSimplex" begin
        for M in (2, 17, 7487)
            for i in (1, 536, typemax(Int32))
                c = rand(Int64)
                d = rand(Float64)
                @test index(DiameterSimplex{M}(d, i, c)) == i
                @test coef(DiameterSimplex{M}(d, i, c)) == mod(c, M)
                @test diam(DiameterSimplex{M}(d, i, c)) == d
            end
        end
        @test_throws DomainError DiameterSimplex{1}(1.0, 1, 1)
        @test_throws DomainError DiameterSimplex{6666}(1.0, 1, 1)
    end

    @testset "index(::Vector), vertices" begin
        st = ReductionState{2}(rand_dist_matrix(10), 5)
        buff = Int[]
        @test vertices(st, Simplex{2}(1, 1), 2) == [3, 2, 1]
        @test vertices(st, Simplex{2}(2, 1), 3) == [5, 3, 2, 1]
        @test vertices(st, Simplex{2}(3, 1), 1) == [3, 2]
        @test vertices(st, Simplex{2}(4, 1), 4) == [6, 5, 4, 2, 1]
        @test vertices(st, Simplex{2}(5, 1), 2) == [5, 2, 1]

        for i in 1:10
            @test index(st, vertices(st, Simplex{2}(i, 1), 5)) == i
        end
    end

    @testset "set_coef" begin
        @test set_coef(Simplex{3}(2, 1), 2) == Simplex{3}(2, 2)
        @test set_coef(Simplex{3}(2, 1), 3) == Simplex{3}(2, 0)
        @test set_coef(DiameterSimplex{5}(1.0, 2, 1), 3) == DiameterSimplex{5}(1.0, 2, 3)
    end

    @testset "arithmetic" begin
        @test Simplex{3}(3, 2) * 2 == Simplex{3}(3, 1)
        @test 2 * Simplex{3}(3, 2) == Simplex{3}(3, 1)
        @test -Simplex{5}(1, 1) == Simplex{5}(1, 4)

        @test inv_mod(Val(2), 1) == 1
        @test_throws DomainError inv_mod(Val(4), 1)
        @test_throws DivideError inv_mod(Val(3), 0)

        for i in 1:16
            @test Simplex{17}(10, i) / i == Simplex{17}(10, 1)
            @test -Simplex{17}(8, i) == Simplex{17}(8, 17 - i)
        end
    end
end
