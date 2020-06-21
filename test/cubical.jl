using Ripserer

using Compat
using StaticArrays

using Ripserer: all_equal_in_dim, boundary, coboundary, face_type, coface_type,
    edges, n_vertices, index

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

    @testset "Vertices, index." begin
        @testset "Basics." begin
            @test vertices(Cubelet{2}(1, rand())) === SVector{4}(4, 3, 2, 1)
            @test vertices(Cubelet{2}(Int128(2), rand())) === SVector{4, Int128}(5, 3, 2, 1)
            @test vertices(Cubelet{1}(3, rand())) === SVector{2}(3, 2)
            @test vertices(Cubelet{0}(Int32(4), rand())) === SVector{1, Int32}(4)
            @test vertices(Cubelet{3}(5, rand())) === SVector{8}(9, 8, 7, 6, 4, 3, 2, 1)
        end

        @testset "Index to vertices and back." begin
            for i in 1:20
                cube = Cubelet{5}(i, rand())
                @test index(vertices(cube)) == i
                @test Cubelet{5}(vertices(cube), diam(cube)) == cube
                @test_throws ArgumentError Cubelet{4}(vertices(cube), rand())
                @test_throws ArgumentError Cubelet{6}(vertices(cube), rand())
            end
        end
    end
    @testset "Printing." begin
        @test sprint(print, Cubelet{1}(1, 1)) == "Cubelet{1}(+[2, 1], 1)"
        @test sprint(print, Cubelet{2}(-1, 1)) == "Cubelet{2}(-[4, 3, 2, 1], 1)"

        @test sprint(Cubelet{2}(1, 1)) do io, sx
            show(io, MIME"text/plain"(), sx)
        end == "2-dim Cubelet(1, 1):\n  +[4, 3, 2, 1]"

        @test sprint(Cubelet{1}(Int128(1), 1)) do io, sx
            show(io, MIME"text/plain"(), sx)
        end == "1-dim Cubelet(1, 1):\n  +Int128[2, 1]"
    end
end

@testset "Cubical" begin
    test_filtration_interface(Cubical, (data1d, data2d, data3d))

    @testset "n_vertices, indices, birth, diam" begin
        for data in (data1d, data2d, data3d)
            filtration = Cubical(data)

            @test n_vertices(filtration) == length(data)
            @test CartesianIndices(filtration) == CartesianIndices(data)
            @test LinearIndices(filtration) == LinearIndices(data)
            @test birth(filtration, 10) == data[10]
            #@test diam(filtration, (10, 9)) == max(data[10], data[9])
        end
    end
end

@testset "Coboundary." begin
    @testset "all_equal_in_dim" begin
        @test all_equal_in_dim(1, [(1, 1), (1, 2), (1, 3), (1, 4)])
        @test !all_equal_in_dim(2, [(1, 1), (1, 1), (1, 2), (1, 2)])
    end

    @testset "Edges have no coboundary in 1d." begin
        cob = Cubelet{2, Float64, Int}[]
        flt = Cubical(data1d)
        cub = simplex(flt, Val(1), (3, 2))
        for c in coboundary(flt, cub)
            push!(cob, c)
        end
        @test isempty(cob)
    end

    @testset "Coboundary of edges in 2d." begin
        cob = Cubelet{2, Int, Int}[]
        flt = Cubical(data2d)
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

    @testset "Diameters of 2-cubelet coboundary in 2×2×2 3d data." begin
        data = zeros(2, 2, 2)
        data[:, 1, 1] .= 1
        data[:, 2, 1] .= 2
        data[:, 1, 2] .= 3
        data[:, 2, 2] .= 4

        @testset "Diameter increasing, cofaces are supersets of original." begin
            cob = Cubelet{2, Float64, Int}[]
            for cube in edges(Cubical(data))
                for c in coboundary(Cubical(data), cube)
                    @test issubset(cube, c)
                    @test diam(cube) ≤ diam(c)
                end
            end
        end
        @testset "Coboundary with all_cofaces=false." begin
            cob = Cubelet{2, Float64, Int}[]
            for cube in edges(Cubical(data))
                for c in coboundary(Cubical(data), cube, Val(false))
                    push!(cob, c)
                end
            end
            @test length(cob) == 6
            @test allunique(cob)
        end
        @testset "All 2-cubelets have the same coboundary." begin
            for cube in [
                Cubelet{2}((1, 2, 3, 4), 2.0),
                Cubelet{2}((1, 2, 5, 6), 4.0),
                Cubelet{2}((1, 3, 5, 7), 3.0),
                Cubelet{2}((2, 4, 6, 8), 4.0),
                Cubelet{2}((3, 4, 7, 8), 4.0),
                Cubelet{2}((5, 6, 7, 8), 4.0),
            ]
                coface = only(coboundary(Cubical(data), cube))
                @test abs(coface) == Cubelet{3}((8, 7, 6, 5, 4, 3, 2, 1), 4.0)
            end
        end
    end
    @testset "Counting cofaces in 3d." begin
        cob = Cubelet{2, Float64, Int}[]
        flt = Cubical(data3d)
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

@testset "Boundary." begin
    @testset "Boundary of an edge in 1d." begin
        bnd = Cubelet{0, Float64, Int}[]
        flt = Cubical(data3d)
        cub = Cubelet{1}((3, 2), 1.0)
        for f in boundary(flt, cub)
            push!(bnd, f)
        end
        @test bnd == [Cubelet{0}((3,), 1.0), -Cubelet{0}((2,), 1.0)]
    end
    @testset "Boundary of 2-cubelet in 2d." begin
        bnd = Cubelet{1, Int, Int}[]
        flt = Cubical(data2d)
        cub = Cubelet{2}((10, 9, 3, 2), 2)
        for f in boundary(flt, cub)
            push!(bnd, f)
        end
        @test bnd == [Cubelet{1}((10, 3), 2), -Cubelet{1}((9, 2), 2),
                      Cubelet{1}((10, 9), 1), -Cubelet{1}((3, 2), 2)]
    end
    @testset "Boundary of 3-cubelet in 3d." begin
        bnd = Cubelet{2, Float64, Int}[]
        flt = Cubical(data3d)
        cub = Cubelet{3}((112, 111, 102, 101, 12, 11, 2, 1), 1.0)
        for f in boundary(flt, cub)
            push!(bnd, f)
        end
        @test bnd == [Cubelet{2}((112, 102, 12, 2), 1.0),
                      -Cubelet{2}((111, 101, 11, 1), 1.0),
                      Cubelet{2}((112, 111, 12, 11), 1.0),
                      -Cubelet{2}((102, 101, 2, 1), 1.0),
                      Cubelet{2}((112, 111, 102, 101), 1.0),
                      -Cubelet{2}((12, 11, 2, 1), 1.0),
        ]
    end
end
