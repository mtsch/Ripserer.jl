using Ripserer: dist_type, dist

"""
    filtration_interface_testset(Filtration, datasets)

Test the filtration interface. Make sure all required methods are defined and that they
return values of types that make sense.
"""
function test_filtration_interface(Filtration, datasets)
    @testset "$(string(Filtration))" begin
        for data in datasets
            flt = Filtration(data)

            T = dist_type(flt)
            S = eltype(flt)

            @test S <: AbstractSimplex

            @test n_vertices(flt) == size(data, 1)
            l, (i, j) = first(edges(flt))
            @test i isa Integer
            @test j isa Integer
            @test l isa T

            @test diam(flt, [3, 2, 1]) isa T
            @test diam(flt, S(one(T), 1, 1), [3, 2], 1) isa T
            @test issparse(flt) isa Bool
            @test issparse(Filtration) isa Bool
            @test issparse(Filtration) == issparse(flt)
        end
    end
end

@testset "rips" begin
    for Filtration in (RipsFiltration, SparseRipsFiltration)
        test_filtration_interface(Filtration, (icosahedron, torus_points(100)))

        @testset "$(string(Filtration))" begin
            @testset "length, dist, threshold, eltype" begin
                flt = Filtration([0 1 2 9;
                                  1 0 3 9;
                                  2 3 0 4;
                                  9 9 4 0], modulus=3)
                @test n_vertices(flt) == 4
                @test dist(flt, 3, 3) == 0
                @test dist(flt, 1, 2) == 1
                @test dist(flt, 1, 3) == 2
                @test dist(flt, 3, 2) == 3
                @test diam(flt, Simplex{3}(3, 3, 1), [3, 2], 1) == 3
                @test flt.threshold == 4
                @test eltype(flt) === Simplex{3, Int}

                flt = Filtration([0 1 2;
                                  1 0 3;
                                  2 3 0], threshold=2, simplex_type=Simplex{5, Int})
                @test n_vertices(flt) == 3
                @test dist(flt, 3, 3) == 0
                @test dist(flt, 1, 2) == 1
                @test dist(flt, 1, 3) == 2
                @test dist(flt, 3, 2) == (issparse(Filtration) ? ∞ : 3)
                @test diam(flt, Simplex{5}(1, 1, 1), [1, 2], 3) == ∞
                @test flt.threshold == 2
                @test eltype(flt) === Simplex{5, Int}
            end
        end
    end
end
