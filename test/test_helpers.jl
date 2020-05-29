module TestHelpers
using Test

using Ripserer
using Ripserer: dist_type, vertex_type, edge_type

"""
    test_indexed_simplex_interface(S, n_vertices)

Test the `IndexedSimplex` interface for a simplex type `S` where a `D`-simplex has
`n_vertices(D)` vertices.
"""
function test_indexed_simplex_interface(S, n_vertices)
    @testset "interface" begin
        for D in (0, 2, 5)
            for T in (Float64, Float32, Int)
                d = rand(T)
                for I in (Int64, Int128)
                    i = rand(I)

                    @test index(S{D}(i, d)) == i
                    @test diam(S{D}(i, d)) == d
                    @test typeof(S{D}(i, d)) ≡ S{D, T, I}
                    @test S{D}(i, d) isa IndexedSimplex{D, T, I}
                    @test S{D}(i, d) isa AbstractSimplex{D, T}
                    D > 0 && @test_throws DomainError S{-D}(i, d)

                    @test -S{D}(i, d) == S{D}(-i, d)
                    @test sign(+S{D}(i, d)) == sign(i)
                    @test sign(-S{D}(i, d)) == -sign(i)

                    @test coface_type(S{D}(i, d)) ≡ S{D + 1, T, I}
                    @test coface_type(typeof(S{D}(i, d))) ≡ S{D + 1, T, I}

                    @test dim(S{D}(i, d)) == D
                    @test abs(S{D}(i, d)) == S{D}(abs(i), d)
                    @test abs(-S{D}(i, d)) == abs(S{D}(i, d))
                end
                # don't want to do this with random int value.
                @test length(vertices(S{D}(10, d))) == n_vertices(D)
                @test length(vertices(S{D}(-10, d))) == n_vertices(D)

                @test begin @inferred vertices(S{D}(10, d)); true end
                @test begin @inferred vertices(S{D}(-10, d)); true end
            end
        end

        for sgn in (1, -1)
            @test isless(S{0}(sgn * 10, 1), S{0}(10, 2))
            @test !isless(S{1}(sgn * 10, 2), S{1}(10, 1))
            @test isless(S{2}(sgn * 11, 1), S{2}(10, 1))
            @test !isless(S{3}(sgn * 9, 1), S{3}(10, 1))
            @test !isless(S{4}(sgn * 10, 1), S{4}(10, 1))
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

            T = dist_type(flt)
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
            @test birth(flt, 1) isa T

            @test diam(flt, (4, 3, 2, 1)) isa Union{T, Missing}

            @test begin @inferred edges(flt); true end
        end
    end
end

export test_simplex_interface, test_filtration_interface
end
