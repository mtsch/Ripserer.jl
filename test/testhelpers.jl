module TestHelpers
using Test

using Ripserer
using Ripserer: vertex_type, edges, edge_type,
    AbstractSimplex, IndexedSimplex, index

"""
    test_indexed_simplex_interface(S, n_vertices)

Test the `IndexedSimplex` interface for a simplex type `S` where a `D`-simplex has
`n_vertices(D)` vertices.
"""
function test_indexed_simplex_interface(S, n_vertices)
    @testset "interface" begin
        for D in (0, 2, 4),
            T in (Float64, Float32, Int),
            I in (Int32, Int64, Int128)

            # Check for overflow.
            try
                vxs_sml = vertices(I(1000), Val(n_vertices(D)))
                vxs_big = vertices(big(1000), Val(n_vertices(D)))
                if vxs_sml ≠ vxs_big
                    continue
                end
            catch InexactError
                continue
            end

            @testset "basics" begin
                d = rand(T)
                i = I(1000)

                @test S{D}(I(i), d) ≡ S{D, T, I}(i, d)
                @test S{D}(I(i), T(1)) ≡ S{D, T, I}(i, 1)

                @test index(S{D}(i, d)) ≡ i
                @test diam(S{D}(i, d)) ≡ d
                @test typeof(S{D}(i, d)) ≡ S{D, T, I}
                @test S{D}(i, d) isa IndexedSimplex{D, T, I}
                D > 0 && @test_throws DomainError S{-D}(i, d)

                @test -S{D}(i, d) == S{D}(-i, d)
                @test sign(+S{D}(i, d)) == sign(i)
                @test sign(-S{D}(i, d)) == -sign(i)

                @test dim(S{D}(i, d)) == D
                @test abs(S{D}(i, d)) == S{D}(abs(i), d)
                @test abs(-S{D}(i, d)) == abs(S{D}(i, d))

                @test eltype(S{D}(i, d)) == I

                for sgn in (1, -1)
                    @test isless(S{D}(sgn * 10, 1), S{D}(10, 2))
                    @test !isless(S{D}(sgn * 10, 2), S{D}(10, 1))
                    @test isless(S{D}(sgn * 11, 1), S{D}(10, 1))
                    @test !isless(S{D}(sgn * 9, 1), S{D}(10, 1))
                    @test !isless(S{D}(sgn * 10, 1), S{D}(10, 1))
                end
            end

            @testset "array interface, vertices" begin
                d = rand(T)
                # don't want to do this with random int value.
                @test eltype(S{D}(I(9), d)) == eltype(vertices(S{D}(I(9), d)))
                @test length(S{D}(I(9), d)) == n_vertices(D)
                @test size(S{D}(-I(9), d)) == (n_vertices(D),)
                @test length(vertices(S{D}(I(10), d))) == n_vertices(D)
                @test length(vertices(S{D}(-I(10), d))) == n_vertices(D)
                @test firstindex(S{D}(I(11), d)) == 1
                @test lastindex(S{D}(I(11), d)) == n_vertices(D)
                @test (1:1000)[S{D}(-I(12), d)] == Int.(vertices(S{D}(I(12), d)))

                @test begin @inferred vertices(S{D}(I(10), d)); true end
                @test begin @inferred vertices(S{D}(-I(10), d)); true end
            end
        end
    end
end

"""
    filtration_interface_testset(Filtration, datasets)

Test the filtration interface. Make sure all required methods are defined and that they
return values of types that make sense.
"""
function test_filtration_interface(Filtration, datasets)
    @testset "AbstractFiltration interface" begin
        for data in datasets
            flt = Filtration(data)

            T = eltype(data)
            V = vertex_type(flt)
            E = edge_type(flt)

            @test V <: AbstractSimplex{0}
            @test V isa DataType
            @test E <: AbstractSimplex{1}
            @test E isa DataType

            @test threshold(flt) isa Union{T, Missing}

            sx = first(edges(flt))
            @test sx isa E
            @test diam(sx) isa T
            @test birth(flt, 1) isa Union{T, Bool}

            #@test diam(flt, (4, 3, 2, 1)) isa Union{T, Missing}

            @test begin @inferred edges(flt); true end
        end
    end
end

export test_indexed_simplex_interface, test_filtration_interface
end
