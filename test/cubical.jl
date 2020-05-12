using Ripserer
using Ripserer: all_equal_in_dim

data1d = cos.(range(0, 4π, length=1000))

data2d = [2 2 2 0 1 1 1 1 1 1 1;
          2 1 2 0 1 2 2 2 2 2 1;
          2 1 2 0 1 2 3 3 3 2 1;
          2 0 2 0 1 2 3 2 3 2 1;
          2 1 2 0 1 2 3 3 3 2 1;
          2 1 2 0 1 2 2 2 2 2 1;
          2 2 2 0 1 1 1 1 1 1 1]

data3d = reshape(fill(1.0, 1000), (10, 10, 10))

using ..TestHelpers: test_indexed_simplex_interface, test_filtration_interface

@testset "Cubelet" begin
    test_indexed_simplex_interface(Cubelet, D -> 2^D)

    @testset "vertices, index" begin
        @test vertices(Cubelet{2}(1, rand())) == (4, 3, 2, 1)
        @test vertices(Cubelet{2}(2, rand())) == (5, 3, 2, 1)
        @test vertices(Cubelet{1}(3, rand())) == (3, 2)
        @test vertices(Cubelet{0}(4, rand())) == (4,)
        @test vertices(Cubelet{3}(5, rand())) == (9, 8, 7, 6, 4, 3, 2, 1)

        for i in 1:20
            sx = Cubelet{5}(i, rand())
            @test index(vertices(sx)) == i
            @test Cubelet{5}(vertices(sx), diam(sx)) == sx
            @test_throws ArgumentError Cubelet{4}(vertices(sx), rand())
            @test_throws ArgumentError Cubelet{6}(vertices(sx), rand())
        end
    end
    @testset "show" begin
        @test sprint(print, Cubelet{1}(1, 1)) == "Cubelet{1}(+(2, 1), 1)"
        @test sprint(print, Cubelet{2}(-1, 1)) == "Cubelet{2}(-(4, 3, 2, 1), 1)"

        @test sprint(Cubelet{2}(1, 1)) do io, sx
            show(io, MIME"text/plain"(), sx)
        end == "2-dim Cubelet(1, 1):\n  +(4, 3, 2, 1)"

        @test sprint(Cubelet{1}(Int128(1), 1)) do io, sx
            show(io, MIME"text/plain"(), sx)
        end == "1-dim Cubelet(1, 1) with Int128 index:\n  +(2, 1)"
    end
end

@testset "CubicalFiltration" begin
    test_filtration_interface(CubicalFiltration, (data1d, data2d, data3d))

    @testset "n_vertices, indices, birth, diam" begin
        for data in (data1d, data2d, data3d)
            flt = CubicalFiltration(data)

            @test n_vertices(flt) == length(data)
            @test CartesianIndices(flt) == CartesianIndices(data)
            @test LinearIndices(flt) == LinearIndices(data)
            @test birth(flt, 10) == data[10]
            @test diam(flt, (10, 9)) == max(data[10], data[9])
        end
    end
end

@testset "coboundary" begin
    @testset "all_equal_in_dim" begin
        @test all_equal_in_dim(1, [(1, 1), (1, 2), (1, 3), (1, 4)])
        @test !all_equal_in_dim(2, [(1, 1), (1, 1), (1, 2), (1, 2)])
    end

    @testset "1d" begin
        cob = Cubelet{2, Float64, Int}[]
        flt = CubicalFiltration(data1d)
        cub = Cubelet{1}((3, 2), diam(flt, (3, 2)))
        for c in coboundary(flt, cub)
            push!(cob, c)
        end
        @test isempty(cob)
    end

    @testset "2d" begin
        cob = Cubelet{2, Int, Int}[]
        flt = CubicalFiltration(data2d)
        cub = Cubelet{1}((10, 9), 1)
        for c in coboundary(flt, cub)
            push!(cob, c)
        end
        @test length(cob) == 2
        @test cob[1] == Cubelet{2}((17, 16, 10, 9), 2)
        @test sign(cob[1]) == 1
        @test cob[2] == -Cubelet{2}((10, 9, 3, 2), 2)
        @test sign(cob[2]) == -1

        cocob = coface_type(eltype(cob))[]
        for c in coboundary(flt, cob[1])
            push!(cocob, c)
        end
        for c in coboundary(flt, cob[2])
            push!(cocob, c)
        end
        @test isempty(cocob)
    end

    @testset "3d" begin
        cob = Cubelet{2, Int, Int}[]
        flt = CubicalFiltration(data3d)
        cub = Cubelet{1}((14, 13), 1)
        for c in coboundary(flt, cub)
            push!(cob, c)
        end
        @test length(cob) == 3

        cocob = coface_type(eltype(cob))[]
        for c in coboundary(flt, cob[1])
            push!(cocob, c)
        end
        @test length(cocob) == 1
    end
end
