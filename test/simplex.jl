using Ripserer
using StaticArrays
using Test

using Ripserer: _binomial, _vertices, coboundary, boundary

struct FakeFiltration<:Ripserer.AbstractFiltration{Int, Int} end
function Ripserer.unsafe_simplex(
    ::Type{Simplex{D, Int, Int}},
    ::FakeFiltration,
    vertices,
    sign,
) where D
    return Simplex{D, Int, Int}(sign * index(vertices), 1)
end
Ripserer.nv(::FakeFiltration) = 20
Ripserer.simplex_type(::Type{FakeFiltration}, D) = Simplex{D, Int, Int}

struct FakeFiltrationWithThreshold<:Ripserer.AbstractFiltration{Int, Int} end
function Ripserer.unsafe_simplex(
    ::Type{Simplex{D, Int, Int}}, ::FakeFiltrationWithThreshold, vertices, sign
) where D
    if maximum(vertices) > 10
        return nothing
    else
        return Simplex{D, Int, Int}(sign * index(vertices), 1)
    end
end
Ripserer.nv(::FakeFiltrationWithThreshold) = 20
Ripserer.simplex_type(::Type{FakeFiltrationWithThreshold}, D) = Simplex{D, Int, Int}

@testset "Internal functions" begin
    @testset "_binomial" begin
        @test all(_binomial(n, Val(k)) == binomial(n, k) for n in 0:1000 for k in 0:7)
    end

    @testset "index" begin
        @test index((3, 2, 1)) === 1
        @test index(Int32.((4, 3, 2, 1))) === Int32(1)
        @test index(SVector(101, 100, 20, 15, 3)) ==
            1 +
            binomial(101 - 1, 5) +
            binomial(100 - 1, 4) +
            binomial(20 - 1, 3) +
            binomial(15 - 1, 2) +
            binomial(3 - 1, 1)
    end

    @testset "_vertices" begin
        for i in 1:3:10, N in 1:5
            @test index(_vertices(i, Val(N))) == i
            @test length(_vertices(i, Val(N))) == N
            @test begin @inferred _vertices(i, Val(N)); true end
        end
    end
end

@testset "Simplex" begin
    @testset "Constructors" begin
        @test Simplex{1}([2, 1], 1) === Simplex{1}(1, 1)
        @test Simplex{1}((1, 3), 0f0) === Simplex{1}(2, 0f0)
        @test Simplex{2}([5, 2, 1], 1.0) === Simplex{2}(5, 1.0)
    end

    @testset "AbstractSimplex interface" begin
        for D in (0, 2, 4),
            T in (Float64, Float32, Int),
            I in (Int32, Int64, Int128)
            @testset "Simplex{$D, $T, $I}" begin
                b = T(9)
                i = I(100)
                sx = Simplex{D}(i, b)

                @testset "Type stuff" begin
                    @test sx ≡ Simplex{D, T, I}(i, b)
                    @test typeof(sx) ≡ Simplex{D, T, I}
                    D > 0 && @test_throws DomainError Simplex{-D}(i, b)
                    @test eltype(sx) == I
                end

                @testset "Getters" begin
                    @test index(sx) == i
                    @test index(-sx) == i
                    @test birth(sx) ≡ b
                end

                @testset "Equality, hashing" begin
                    @test sx == sx
                    @test sx != Simplex{D + 1}(i, b)
                    @test sx == Simplex{D}(i, b + 1)
                    @test isequal(sx, Simplex{D}(i, b + 1))
                    @test hash(sx) == hash(Simplex{D}(i, b + 1))
                    @test hash(sx) == hash(index(sx))
                end

                @testset "Signs" begin
                    @test +sx === sx
                    @test sx == -sx
                    @test sx ≡ -(-sx)
                    @test sign(sx) == 1
                    @test sign(-sx) == -1
                    @test abs(sx) ≡ sx
                    @test abs(-sx) ≡ sx
                    @test dim(sx) == D
                end

                @testset "Ordering" begin
                    @test sx < Simplex{D}(I(i + 1), b + 1)
                    @test sx > Simplex{D}(I(i - 1), b - 1)

                    @test sx < Simplex{D}(I(i - 1), b)
                    @test sx > Simplex{D}(I(i + 1), b)
                end

                @testset "Array interface, vertices" begin
                    verts = vertices(sx)

                    @test eltype(sx) == eltype(verts)
                    @test length(sx) == length(verts) == D + 1
                    @test size(sx) == (D + 1,)
                    @test firstindex(sx) == 1
                    @test lastindex(sx) == D + 1

                    @test Simplex{D}(verts, birth(sx)) ≡ sx

                    for (i, v) in enumerate(sx)
                        @test v == verts[i]
                    end

                    @test begin @inferred vertices(sx); true end
                end

                @testset "Printing" begin
                    @test sprint(show, sx) ==
                        "+Simplex{$D}($(vertices(sx)), $(b))"
                    @test sprint(show, -sx) ==
                        "-Simplex{$D}($(vertices(sx)), $(b))"
                    @test sprint((i, s) -> show(i, MIME"text/plain"(), s), sx) ==
                        "$D-dimensional Simplex(index=$i, birth=$b):\n  +$(vertices(sx))"
                end
            end
        end
    end

    @testset "Coboundary" begin
        @testset "number of cofacets" begin
            for dim in 1:5
                simplex = Simplex{dim}(10, 1)
                simplex_vxs = vertices(simplex)
                cob = []
                for cofacet in coboundary(FakeFiltration(), simplex)
                    push!(cob, cofacet)
                end
                @test length(cob) == 20 - dim - 1
                @test all(isequal(1), birth.(cob))
                @test all(issubset(simplex_vxs, vertices(c)) for c in cob)
                @test issorted(cob, by=index, rev=true)
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
            for dim in 1:5
                simplex = Simplex{dim}(10, 1)
                simplex_vxs = vertices(simplex)
                cob = []
                for cofacet in coboundary(FakeFiltrationWithThreshold(), simplex)
                    push!(cob, cofacet)
                end
                @test length(cob) == 10 - dim - 1
                @test all(isequal(1), birth.(cob))
                @test all(issubset(simplex_vxs, vertices(c)) for c in cob)
                @test issorted(cob, by=index, rev=true)
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

    @testset "Boundary" begin
        @testset "number of facets" begin
            for dim in 1:5
                simplex = Simplex{dim}(10, 1)
                simplex_vxs = vertices(simplex)
                bnd = []
                for facet in boundary(FakeFiltration(), simplex)
                    push!(bnd, facet)
                end
                @test length(bnd) == dim + 1
                @test all(isequal(1), birth.(bnd))
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
