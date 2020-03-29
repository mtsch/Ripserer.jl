using Ripserer: Simplex, DiameterSimplex, index, coef, diam, get_vertices!

@testset "simplices" begin
    @testset "Simplex" begin
        for modulus in (2, 17, 7487)
            for i in (1, 536, typemax(Int32))
                c = rand(Int64)
                @test Simplex(i, c, modulus) === Simplex{modulus}(i, c)
                @test index(Simplex(i, c, modulus)) == i
                @test coef(Simplex(i, c, modulus)) == mod(c, modulus)
            end
        end
        @test_throws DomainError Simplex(1, 1, 4)
        @test_throws DomainError Simplex(1, 1, 7497)
    end

    @testset "DiameterSimplex" begin
        for modulus in (2, 17, 7487)
            for i in (1, 536, typemax(Int32))
                c = rand(Int64)
                d = rand(Float64)
                @test DiameterSimplex(d, i, c, modulus) === DiameterSimplex{modulus}(d, i, c)
                @test index(DiameterSimplex(d, i, c, modulus)) == i
                @test coef(DiameterSimplex(d, i, c, modulus)) == mod(c, modulus)
                @test diam(DiameterSimplex(d, i, c, modulus)) == d
            end
        end
        @test_throws DomainError DiameterSimplex(1.0, 1, 1, 1)
        @test_throws DomainError DiameterSimplex(1.0, 1, 1, 6666)
    end

    @testset "index, get_vertices!" begin
        buff = Int[]
        @test get_vertices!(buff, Simplex(1, 1, 2), 2, 10, binomial) == [3, 2, 1]
        @test get_vertices!(buff, Simplex(2, 1, 2), 2, 10, binomial) == [4, 2, 1]
        @test get_vertices!(buff, Simplex(3, 1, 2), 2, 10, binomial) == [4, 3, 1]
        @test get_vertices!(buff, Simplex(4, 1, 2), 2, 10, binomial) == [4, 3, 2]
        @test get_vertices!(buff, Simplex(5, 1, 2), 2, 10, binomial) == [5, 2, 1]

        for i in 1:10
            @test index(get_vertices!(buff, Simplex(i, 1, 2), 0, 10, binomial), binomial) == i
        end
    end
end
