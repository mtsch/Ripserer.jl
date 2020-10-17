using Compat
using Random
using Ripserer
using StaticArrays
using Test
using TupleTools

using Ripserer:
    _one_hot,
    _cubemap,
    _from_cubemap,
    _to_cubemap,
    nv,
    coboundary,
    boundary,
    edges,
    births,
    chain_element_type,
    CubicalChainElement

@testset "CubeMap" begin
    @testset "cubemap" begin
        @testset "1d" begin
            data = Float64[9, 1, 8, 2, 1, 3]
            @test _cubemap(data) == [9, 9, 1, 8, 8, 8, 2, 2, 1, 3, 3]
        end
        @testset "2d" begin
            data = [1 2; 3 4]
            @test _cubemap(data) == [1 2 2; 3 4 4; 3 4 4]

            data = [1; 2; 3]
            @test _cubemap(data) == [1; 2; 2; 3; 3]
        end
        @testset "nd" begin
            for D in 3:6
                data = ones(ntuple(identity, D)...)
                map = _cubemap(data)

                @test map == ones(ntuple(i -> 2i - 1, D)...)
            end
        end
    end

    @testset "from/to_cubemap" begin
        @test _from_cubemap(CartesianIndex(2, 2, 2), Val(8)) == SVector(
            CartesianIndex(1, 1, 1),
            CartesianIndex(2, 1, 1),
            CartesianIndex(1, 2, 1),
            CartesianIndex(2, 2, 1),
            CartesianIndex(1, 1, 2),
            CartesianIndex(2, 1, 2),
            CartesianIndex(1, 2, 2),
            CartesianIndex(2, 2, 2),
        )

        for vertices in (
            (CartesianIndex(5, 4),),
            (CartesianIndex(1024, 1072), CartesianIndex(1025, 1072)),
            (
                CartesianIndex(1024, 1072, 15),
                CartesianIndex(1025, 1072, 15),
                CartesianIndex(1024, 1073, 15),
                CartesianIndex(1025, 1073, 15),
                CartesianIndex(1024, 1072, 14),
                CartesianIndex(1025, 1072, 14),
                CartesianIndex(1024, 1073, 14),
                CartesianIndex(1025, 1073, 14),
            ),
        )
            @test _from_cubemap(_to_cubemap(vertices), Val(length(vertices))) ==
                  sort(SVector(vertices))
            root = _to_cubemap(vertices)
            @test begin
                _from_cubemap(root, Val(length(vertices)))
                true
            end
        end

        @test_throws ArgumentError _from_cubemap(CartesianIndex(2, 2, 1), Val(8))
    end
end

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
        for (D, K) in ((0, 5), (1, 1), (2, 3), (4, 4)), T in (Float32, Int)
            @testset "Cube{$D, $T, $K}" begin
                d = rand(T)
                root = CartesianIndex{K}(ntuple(_ -> 2rand(1:100) - 1, Val(K)))
                dirs = Tuple(sort!(shuffle(1:K)[1:D]))
                for i in dirs
                    root += _one_hot(i, Val(K))
                end
                c = Cube{D}(root, d)

                @testset "Type stuff" begin
                    @test c ≡ Cube{D,T,K}(root, d)
                    @test typeof(c) ≡ Cube{D,T,K}
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

                    @test c < Cube{D}(root - 2_one_hot(rand(1:K), Val(K)), d)
                    @test c > Cube{D}(root + 2_one_hot(rand(1:K), Val(K)), d)
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

                    @test begin
                        @inferred vertices(c)
                        true
                    end
                end
            end
        end
    end

    @testset "Vertices" begin
        @test vertices(Cube{0}(((1, 1),), 1)) == [CartesianIndex(1, 1)]
        @test vertices(Cube{1}(((1, 1), (1, 2)), 1)) == CartesianIndex.([(1, 1), (1, 2)])
        @test vertices(Cube{4}(CartesianIndex(10, 14, 15, 22, 2), 1)) ==
              CartesianIndex.([
            (5, 7, 8, 11, 1),
            (6, 7, 8, 11, 1),
            (5, 8, 8, 11, 1),
            (6, 8, 8, 11, 1),
            (5, 7, 8, 12, 1),
            (6, 7, 8, 12, 1),
            (5, 8, 8, 12, 1),
            (6, 8, 8, 12, 1),
            (5, 7, 8, 11, 2),
            (6, 7, 8, 11, 2),
            (5, 8, 8, 11, 2),
            (6, 8, 8, 11, 2),
            (5, 7, 8, 12, 2),
            (6, 7, 8, 12, 2),
            (5, 8, 8, 12, 2),
            (6, 8, 8, 12, 2),
        ])
    end
end

@testset "Cubical" begin
    for K in 1:4
        data = rand(Float64, ntuple(x -> 10 + x, Val(K)))
        filtration = Cubical(data)

        @test nv(filtration) == length(data)
        @test births(filtration) == data
        @test size(filtration.cubemap) == size(data) .* 2 .- 1
        @test vertices(filtration) == CartesianIndices(data)
    end
end

@testset "Coboundary" begin
    @testset "Edges have no coboundary in 1d" begin
        data = cos.(range(0, 4π; length=1000))

        cob = Cube{2,Float64,1}[]
        flt = Cubical(data)
        cub = simplex(flt, Val(1), (CartesianIndex(3), CartesianIndex(2)))
        for c in coboundary(flt, cub)
            push!(cob, c)
        end
        @test isempty(cob)
    end

    @testset "Coboundary of edges in 2d" begin
        data = [
            2 2 2 0 1 1 1 1 1 1 1
            2 1 2 0 1 2 2 2 2 2 1
            2 1 2 0 1 2 3 3 3 2 1
            2 0 2 0 1 2 3 2 3 2 1
            2 1 2 0 1 2 3 3 3 2 1
            2 1 2 0 1 2 2 2 2 2 1
            2 2 2 0 1 1 1 1 1 1 1
        ]
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
        @test issorted(cob; by=index, rev=true)

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
                @test issorted(cob; by=index, rev=true)
                @test all(c -> issubset(cube, c), cob)
                @test all(c -> birth(cube) ≤ birth(c), cob)
            end
        end
        @testset "All 2-cubes have the same coboundary" begin
            for cube in [
                Cube{2}(((1, 1, 1), (2, 1, 1), (1, 2, 1), (2, 2, 1)), 2.0),
                Cube{2}(((1, 1, 1), (2, 1, 1), (1, 1, 2), (2, 1, 2)), 4.0),
                Cube{2}(((1, 1, 1), (1, 2, 1), (1, 1, 2), (1, 2, 2)), 3.0),
                Cube{2}(((2, 1, 1), (2, 2, 1), (2, 1, 2), (2, 2, 2)), 4.0),
                Cube{2}(((1, 2, 1), (2, 2, 1), (1, 2, 2), (2, 2, 2)), 4.0),
                Cube{2}(((1, 1, 2), (2, 1, 2), (1, 2, 2), (2, 2, 2)), 4.0),
            ]
                cofacet = only(coboundary(Cubical(data), cube))
                @test abs(cofacet) == Cube{3}(
                    (
                        (2, 2, 2),
                        (1, 2, 2),
                        (2, 1, 2),
                        (1, 1, 2),
                        (2, 2, 1),
                        (1, 2, 1),
                        (2, 1, 1),
                        (1, 1, 1),
                    ),
                    4.0,
                )
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
        @test issorted(cob; by=index, rev=true)

        cocob = []
        for c in coboundary(flt, cob[1])
            push!(cocob, c)
        end
        @test length(cocob) == 1
    end
end

@testset "Boundary" begin
    data = [
        1 2
        3 4
    ]
    cube = Cube{2}(
        (
            CartesianIndex(1, 1),
            CartesianIndex(2, 1),
            CartesianIndex(1, 2),
            CartesianIndex(2, 2),
        ),
        4,
    )

    bnd = []
    for f in boundary(Cubical(data), cube)
        push!(bnd, f)
    end
    @test length(bnd) == 4
    @test issorted(bnd; by=index)
    @test sort(birth.(bnd)) == [2, 3, 4, 4]
end

@testset "CubicalChainElement" begin
    for C in (Cube{2,Float64,3}, Cube{3,Int,5})
        @test @inferred(chain_element_type(C, Mod{2})) == CubicalChainElement{C}
        @test_throws ErrorException chain_element_type(C, Mod{251})
        @test_throws ErrorException chain_element_type(C, Mod{257})
        @test_throws ErrorException chain_element_type(C, Rational{Int})
    end
end

@testset "ripserer" begin
    @testset "1D" begin
        data = [1, 0, 1, 2, 3, 4, 3, 2, 3, 2, 1, 2]
        d0, d1, d2 = ripserer(Cubical(data); dim_max=2)

        @test d0 == [(2, 3), (1, 4), (0, Inf)]
        @test d1 == []
        @test d2 == []
    end
    @testset "1D representatives" begin
        n = 1000
        x = range(0, 1; length=n)
        curve = sin.(2π * 5x) .* x

        d0, _ = ripserer(Cubical(curve); reps=true)

        for int in d0
            birth_sx = birth_simplex(int)
            @test curve[only(birth_sx)] == birth(int) == birth(birth_sx)
        end
        @test sort!(only.(birth_simplex.(d0))) ==
              CartesianIndex.([1, 157, 354, 552, 752, 951])
    end
    @testset "2D image" begin
        data = [
            0 0 0 0 0
            0 2 2 2 0
            0 2 1 2 0
            0 2 2 2 0
            0 0 0 0 0
        ]

        d0, d1, d2 = ripserer(Cubical(data); reps=true, dim_max=2)

        @test d0 == [(1, 2), (0, Inf)]
        @test d1 == [(0, 2)]
        @test d2 == []

        @test vertices(only(representative(d0[1]))) == SVector(CartesianIndex(3, 3))
        @test sort(vertices.(representative(d0[2]))) ==
              sort(SVector.(vec(CartesianIndices(data))))
    end
    @testset "3D image" begin
        # Cube with hole in the middle.
        data = zeros(5, 5, 5)
        data[2, 2:4, 2:4] .= 1
        data[3, :, :] .= [0 0 0 0 0; 0 1 1 1 0; 0 1 0 1 0; 0 1 1 1 0; 0 0 0 0 0]
        data[4, 2:4, 2:4] .= 1

        d0, d1, d2 = ripserer(Cubical(data); dim_max=2)

        @test d0 == [(0, 1), (0, Inf)]
        @test d1 == []
        @test d2 == [(0, 1)]
    end
    @testset "Thresholding" begin
        data = [
            1 1 1
            1 2 1
            1 1 1
        ]
        d0, d1 = ripserer(Cubical(data; threshold=1))
        @test d0 == [(1, Inf)]
        @test d1 == [(1, Inf)]
    end
    @testset "Homology and explicit cohomology" begin
        data = zeros(5, 5, 5)
        data[2, 2:4, 2:4] .= 1
        data[3, :, :] .= [0 0 0 0 0; 0 1 1 1 0; 0 1 0 1 0; 0 1 1 1 0; 0 0 0 0 0]
        data[4, 2:4, 2:4] .= 1

        c = Cubical(data)
        _, hom_imp1, hom_imp2 = ripserer(c; alg=:homology, implicit=true, dim_max=2)
        _, hom_exp1, hom_exp2 = ripserer(c; alg=:homology, implicit=false, dim_max=2)
        _, hom_ass1, hom_ass2 = ripserer(c; alg=:involuted, implicit=false, dim_max=2)
        _, coh_imp1, coh_imp2 = ripserer(c; implicit=true, reps=true, dim_max=2)
        _, coh_exp1, coh_exp2 = ripserer(c; implicit=true, reps=true, dim_max=2)

        @test hom_imp1 == hom_exp1 == hom_ass1 == coh_imp1 == coh_exp1
        @test hom_imp2 == hom_exp2 == hom_ass2 == coh_imp2 == coh_exp2
        @test representative.(hom_imp1) == representative.(hom_exp1)
        @test representative.(hom_imp1) == representative.(hom_ass1)
        @test representative.(hom_imp2) == representative.(hom_exp2)
        @test representative.(hom_imp2) == representative.(hom_ass2)
        @test representative.(coh_imp1) == representative.(coh_exp1)
        @test representative.(coh_imp2) == representative.(coh_exp2)
    end
end
