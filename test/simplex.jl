using Ripserer
using StaticArrays

using Ripserer: small_binomial, boundary, coboundary, index
using ..TestHelpers: test_indexed_simplex_interface

struct FakeFiltration<:Ripserer.AbstractFiltration end
function Ripserer.unsafe_simplex(::FakeFiltration, ::Val{D}, vertices, sign=1) where D
    return Simplex{D, Int, Int}(sign * index(vertices), 1)
end
Ripserer.n_vertices(::FakeFiltration) = 20
Ripserer.simplex_type(::FakeFiltration, D) = Simplex{D, Int, Int}

struct FakeFiltrationWithThreshold<:Ripserer.AbstractFiltration end
function Ripserer.unsafe_simplex(
    ::FakeFiltrationWithThreshold, ::Val{D}, vertices, sign=1
) where D
    if maximum(vertices) > 10
        return nothing
    else
        return Simplex{D, Int, Int}(sign * index(vertices), 1)
    end
end
Ripserer.n_vertices(::FakeFiltrationWithThreshold) = 20
Ripserer.simplex_type(::FakeFiltrationWithThreshold, D) = Simplex{D, Int, Int}

@testset "Binomials" begin
    @test all(binomial(n, k) == small_binomial(n, Val(k)) for n in 0:1000 for k in 0:7)
end

@testset "Simplex" begin
    test_indexed_simplex_interface(Simplex, D->D+1)

    @testset "vertices, index" begin
        @test vertices(Simplex{2}(1, rand())) ≡ SVector{3}(3, 2, 1)
        @test vertices(Simplex{3}(Int32(2), rand())) ≡ SVector{4, Int32}(5, 3, 2, 1)
        @test vertices(Simplex{1}(Int128(3), rand())) ≡ SVector{2, Int128}(3, 2)
        @test vertices(Simplex{4}(4, rand())) ≡ SVector{5}(6, 5, 4, 2, 1)
        @test vertices(Simplex{2}(5, rand())) ≡ SVector{3}(5, 2, 1)

        for i in 1:2:20
            for I in (Int64, Int32, Int128)
                sx = Simplex{5}(I(i), rand())
                @test index(vertices(sx)) ≡ I(i)
                @test Simplex{5}(vertices(sx), diam(sx)) ≡ sx
                @test_throws ArgumentError Simplex{4}(vertices(sx), rand())
                @test_throws ArgumentError Simplex{6}(vertices(sx), rand())
            end
        end
    end
    @testset "show" begin
        @test sprint(print, Simplex{1}(1, 1)) == "Simplex{1}(+[2, 1], 1)"
        @test sprint(print, Simplex{2}(-1, 1)) == "Simplex{2}(-[3, 2, 1], 1)"

        @test sprint(Simplex{2}(1, 1)) do io, sx
            show(io, MIME"text/plain"(), sx)
        end == "2-dim Simplex(1, 1):\n  +[3, 2, 1]"

        @test sprint(Simplex{1}(Int128(1), 1)) do io, sx
            show(io, MIME"text/plain"(), sx)
        end == "1-dim Simplex(1, 1):\n  +Int128[2, 1]"

    end
    @testset "coboundary" begin
        @testset "number of cofacets" begin
            for dim in 1:10
                simplex = Simplex{dim}(10, 1)
                simplex_vxs = vertices(simplex)
                cob = []
                for cofacet in coboundary(FakeFiltration(), simplex)
                    push!(cob, cofacet)
                end
                @test length(cob) == 20 - dim - 1
                @test all(isequal(1), diam.(cob))
                @test all(issubset(simplex_vxs, vertices(c)) for c in cob)
            end
        end
        @testset "number of cofacets, all_cofacets=false" begin
            for dim in 1:5
                cob = Simplex{dim+1, Int, Int}[]
                for idx in binomial(20, dim+1):-1:1
                    simplex = Simplex{dim}(idx, 1)
                    for cofacet in coboundary(FakeFiltration(), simplex, Val(false))
                        push!(cob, cofacet)
                    end
                end
                @test sort(index.(abs.(cob))) == 1:binomial(20, dim+2)
            end
        end
        @testset "number of cofacets, thresholding" begin
            for dim in 1:8 # at 9, simplex with index 10 becomes invalid.
                simplex = Simplex{dim}(10, 1)
                simplex_vxs = vertices(simplex)
                cob = []
                for cofacet in coboundary(FakeFiltrationWithThreshold(), simplex)
                    push!(cob, cofacet)
                end
                @test length(cob) == 10 - dim - 1
                @test all(isequal(1), diam.(cob))
                @test all(issubset(simplex_vxs, vertices(c)) for c in cob)
            end
        end
        @testset "type stability" begin
            for dim in 1:10
                simplex = Simplex{dim}(10, 1)
                @test begin @inferred coboundary(FakeFiltration(), simplex); true end
                cobiter = coboundary(FakeFiltration(), simplex)
                @static if VERSION ≥ v"1.1.0"
                    @test begin @inferred Union{
                        Nothing,
                        Tuple{Simplex{dim + 1, Int, Int}, Tuple{Int, Int}},
                    } iterate(cobiter); true end
                    @test begin @inferred Union{
                        Nothing,
                        Tuple{Simplex{dim + 1, Int, Int}, Tuple{Int, Int}},
                    } iterate(cobiter, (13, dim+1)); true end
                end
            end
        end
    end

    @testset "boundary" begin
        @testset "number of facets" begin
            for dim in 1:10
                simplex = Simplex{dim}(10, 1)
                simplex_vxs = vertices(simplex)
                bnd = []
                for facet in boundary(FakeFiltration(), simplex)
                    push!(bnd, facet)
                end
                @test length(bnd) == dim + 1
                @test all(isequal(1), diam.(bnd))
                @test all(issubset(vertices(f), simplex_vxs) for f in bnd)
            end
        end
        @testset "type stability" begin
            for dim in 1:10
                simplex = Simplex{dim}(10, 1)
                @test begin @inferred boundary(FakeFiltration(), simplex); true end
                biter = boundary(FakeFiltration(), simplex)
                @static if VERSION ≥ v"1.1.0"
                    @test begin @inferred Union{
                        Nothing,
                        Tuple{Simplex{dim - 1, Int, Int}, Int},
                    } iterate(biter); true end
                    @test begin @inferred Union{
                        Nothing,
                        Tuple{Simplex{dim - 1, Int, Int}, Int},
                    } iterate(biter, 1); true end
                end
            end
        end
    end
end
