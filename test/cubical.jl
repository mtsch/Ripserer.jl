using Compat
using Random
using Ripserer
using Test

using Ripserer: one_hot, from_cubemap, to_cubemap, n_vertices,
    coboundary, boundary, edges,
    chain_element_type, CubicalChainElement

@testset "Cube" begin
    @testset "Constructors" begin
        @test Cube{1}([CartesianIndex(1, 1), CartesianIndex(1, 2)], 1) ==
            Cube{1}(CartesianIndex(1, 2), 1)
        @test Cube{1}((CartesianIndex(2, 1), CartesianIndex(1, 1)), 1) ==
            Cube{1}(CartesianIndex(2, 1), 1)
        @test Cube{2}([(1, 1), (1, 2), (2, 2), (2, 1)], 1) ==
            Cube{2}(CartesianIndex(2, 2), 1)
    end

    @testset "Randomized tests for AbstractSimplex interface" begin
        for (D, K) in ((0, 5), (1, 1), (2, 3), (4, 4)), T in (Float64, Float32, Int)
            @testset "Cube{$D, $T, $K}" begin
                d = rand(T)
                root = CartesianIndex{K}(ntuple(_ -> 2rand(1:100) - 1, Val(K)))
                dirs = Tuple(sort!(shuffle(1:K)[1:D]))
                for i in dirs
                    root += one_hot(i, Val(K))
                end
                c = Cube{D}(root, d)

                @testset "Type stuff" begin
                    @test c ≡ Cube{D, T, K}(root, d)
                    @test typeof(c) ≡ Cube{D, T, K}
                    D > 0 && @test_throws DomainError Cube{-D}(root, d)
                    @test eltype(c) == CartesianIndex{K}
                end

                @testset "Signs, getters and equality" begin
                    @test index(c) == root
                    @test c == c
                    @test birth(c) ≡ d
                    @test sign(c) == 1
                    @test -c == c
                    @test abs(c) == c
                    @test dim(c) == D
                end

                @testset "Ordering" begin
                    root2 = CartesianIndex{K}(ntuple(_ -> 2rand(1:100) - 1, Val(K)))
                    @test c < Cube{D}(root2, d + 1)
                    @test c > Cube{D}(root2, d - 1)

                    @test c < Cube{D}(root - 2one_hot(rand(1:K), Val(K)), d)
                    @test c > Cube{D}(root + 2one_hot(rand(1:K), Val(K)), d)
                end

                @testset "Array interface, vertices" begin
                    verts = vertices(c)

                    @test eltype(c) == eltype(verts)
                    @test length(c) == length(verts) == 2^D
                    @test size(c) == (2^D,)
                    @test firstindex(c) == 1
                    @test lastindex(c) == 2^D

                    @test Cube{D}(verts, birth(c)) == c

                    for (i, v) in enumerate(c)
                        @test v == verts[i]
                    end

                    @test begin @inferred vertices(c); true end
                end
            end
        end
    end

    @testset "Vertices" begin
        @test vertices(Cube{0}(((1, 1),), 1)) == [CartesianIndex(1, 1)]
        @test vertices(Cube{1}(((1, 1), (1, 2)), 1)) == CartesianIndex.([(1, 1), (1, 2)])
        @test vertices(Cube{4}(CartesianIndex(10, 14, 15, 22, 2), 1)) ==
            CartesianIndex.([
                (5, 7, 8, 11, 1), (6, 7, 8, 11, 1), (5, 8, 8, 11, 1), (6, 8, 8, 11, 1),
                (5, 7, 8, 12, 1), (6, 7, 8, 12, 1), (5, 8, 8, 12, 1), (6, 8, 8, 12, 1),
                (5, 7, 8, 11, 2), (6, 7, 8, 11, 2), (5, 8, 8, 11, 2), (6, 8, 8, 11, 2),
                (5, 7, 8, 12, 2), (6, 7, 8, 12, 2), (5, 8, 8, 12, 2), (6, 8, 8, 12, 2),
            ])
    end
end

@testset "Cubical" begin
    for K in 1:4
        data = rand(Float64, ntuple(x -> 10 + x, Val(K)))
        filtration = Cubical(data)

        @test n_vertices(filtration) == length(data)
        @test birth(filtration) == data
        @test size(filtration.cubemap) == size(data) .* 2 .- 1
        @test vertices(filtration) == CartesianIndices(data)
    end
end

@testset "Coboundary" begin
    @testset "Edges have no coboundary in 1d" begin
        data = cos.(range(0, 4π, length=1000))

        cob = Cube{2, Float64, 1}[]
        flt = Cubical(data)
        cub = simplex(flt, Val(1), (CartesianIndex(3), CartesianIndex(2)))
        for c in coboundary(flt, cub)
            push!(cob, c)
        end
        @test isempty(cob)
    end

    @testset "Coboundary of edges in 2d" begin
        data = [2 2 2 0 1 1 1 1 1 1 1;
                2 1 2 0 1 2 2 2 2 2 1;
                2 1 2 0 1 2 3 3 3 2 1;
                2 0 2 0 1 2 3 2 3 2 1;
                2 1 2 0 1 2 3 3 3 2 1;
                2 1 2 0 1 2 2 2 2 2 1;
                2 2 2 0 1 1 1 1 1 1 1]
        cob = []
        flt = Cubical(data)
        cub = Cube{1}([(3, 2), (2, 2)], 1)
        for c in coboundary(flt, cub)
            push!(cob, c)
        end
        @test length(cob) == 2
        @test cob[1] == Cube{2}([(3, 3), (2, 3), (3, 2), (2, 2)], 2)
        @test sign(cob[1]) == 1
        @test cob[2] == Cube{2}(((3, 1), (2, 1), (3, 2), (2, 2)), 2)
        @test sign(cob[2]) == 1
        @test issorted(cob, by=index, rev=true)

        cocob = []
        for c in coboundary(flt, cob[1])
            push!(cocob, c)
        end
        for c in coboundary(flt, cob[2])
            push!(cocob, c)
        end
        @test isempty(cocob)
    end

    @testset "births of 2-cubelet coboundary in 2×2×2 3d data" begin
        data = zeros(2, 2, 2)
        data[:, 1, 1] .= 1
        data[:, 2, 1] .= 2
        data[:, 1, 2] .= 3
        data[:, 2, 2] .= 4

        @testset "Birth increasing, cofacets are supersets of original" begin
            for cube in edges(Cubical(data))
                cob = []
                for c in coboundary(Cubical(data), cube)
                    push!(cob, c)
                end
                @test issorted(cob, by=index, rev=true)
                @test all(c -> issubset(cube, c), cob)
                @test all(c -> birth(cube) ≤ birth(c), cob)
            end
        end
        @testset "All 2-cubes have the same coboundary" begin
            for cube in [
                Cube{2}(((1,1,1), (2,1,1), (1,2,1), (2,2,1)), 2.0),
                Cube{2}(((1,1,1), (2,1,1), (1,1,2), (2,1,2)), 4.0),
                Cube{2}(((1,1,1), (1,2,1), (1,1,2), (1,2,2)), 3.0),
                Cube{2}(((2,1,1), (2,2,1), (2,1,2), (2,2,2)), 4.0),
                Cube{2}(((1,2,1), (2,2,1), (1,2,2), (2,2,2)), 4.0),
                Cube{2}(((1,1,2), (2,1,2), (1,2,2), (2,2,2)), 4.0),
            ]
                cofacet = only(coboundary(Cubical(data), cube))
                @test abs(cofacet) == Cube{3}((
                    (2,2,2), (1,2,2), (2,1,2), (1,1,2), (2,2,1), (1,2,1), (2,1,1), (1,1,1)
                ), 4.0)
            end
        end
    end
    @testset "Counting cofacets in 3d" begin
        data = reshape(fill(1.0, 1000), (10, 10, 10))
        cob = []
        flt = Cubical(data)
        cub = Cube{1}(((1, 5, 1), (1, 5, 2)), 1)
        for c in coboundary(flt, cub)
            push!(cob, c)
        end
        @test length(cob) == 3
        @test issorted(cob, by=index, rev=true)

        cocob = []
        for c in coboundary(flt, cob[1])
            push!(cocob, c)
        end
        @test length(cocob) == 1
    end
end

@testset "Boundary" begin
    data = [1 2;
            3 4]
    cube = Cube{2}((
        CartesianIndex(1, 1),
        CartesianIndex(2, 1),
        CartesianIndex(1, 2),
        CartesianIndex(2, 2),
    ), 4)

    bnd = []
    for f in boundary(Cubical(data), cube)
        push!(bnd, f)
    end
    @test length(bnd) == 4
    @test sort(birth.(bnd)) == [2, 3, 4, 4]
end

@testset "CubicalChainElement" begin
    for C in (Cube{2, Float64, 3}, Cube{3, Int, 5})
        @test @inferred(chain_element_type(C, Mod{2})) == CubicalChainElement{C}
        @test_throws ErrorException chain_element_type(C, Mod{251})
        @test_throws ErrorException chain_element_type(C, Mod{257})
        @test_throws ErrorException chain_element_type(C, Rational{Int})
    end
end
