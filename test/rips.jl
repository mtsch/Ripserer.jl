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
            S = edge_type(flt)

            @test S <: AbstractSimplex
            @test S isa DataType

            @test n_vertices(flt) == size(data, 1)
            sx = first(edges(flt))
            @test sx isa S
            @test diam(sx) isa T

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
        test_filtration_interface(Filtration, (icosahedron, torus(100)))

        @testset "$(string(Filtration))" begin
            @testset "no threshold, modulus=3" begin
                flt = Filtration([0 1 2 9;
                                  1 0 3 9;
                                  2 3 0 4;
                                  9 9 4 0], modulus=3)
                @test n_vertices(flt) == 4
                @test dist(flt, 3, 3) == 0
                @test dist(flt, 1, 2) == 1
                @test dist(flt, 1, 3) == 2
                @test dist(flt, 3, 2) == 3
                @test diam(flt, Simplex{1, 3}(3, 3, 1), [3, 2], 1) == 3
                @test threshold(flt) == (issparse(Filtration) ? ∞ : 4)
                @test edge_type(flt) === Simplex{1, 3, Int, UInt}
                @test birth(flt, 1) == 0
                @test birth(flt, 4) == 0
            end
            @testset "threshold, edge_type" begin
                flt = Filtration([1 1 2;
                                  1 2 3;
                                  2 3 3];
                                 threshold=2,
                                 edge_type=Simplex{1, 5, Int, UInt})

                @test n_vertices(flt) == 3
                @test dist(flt, 3, 3) == 0
                @test dist(flt, 1, 2) == 1
                @test dist(flt, 1, 3) == 2
                @test dist(flt, 3, 2) == (issparse(Filtration) ? ∞ : 3)
                @test diam(flt, Simplex{2, 5}(1, 1, 1), [1, 2], 3) == ∞
                @test threshold(flt) == (issparse(Filtration) ? ∞ : 2)
                @test edge_type(flt) === Simplex{1, 5, Int, UInt}

                @test birth(flt, 1) == 1
                @test birth(flt, 2) == 2
                @test birth(flt, 3) == 3
            end
            @testset "errors" begin
                @test_throws ArgumentError Filtration([1 2 3;
                                                       4 5 6;
                                                       7 8 9])
                dists = [0 1 2; 1 0 3; 2 3 0]
                @test_throws ArgumentError Filtration(dists, modulus=4)
                @test_throws ArgumentError Filtration(dists, modulus=-2)
                @test_throws ArgumentError Filtration(dists, edge_type=Int)
                @test_throws ArgumentError Filtration(dists, edge_type=Simplex{2,2,Int,UInt})
                @test_throws TypeError Filtration(dists, edge_type=Simplex{1, 2, Int})
            end
        end
    end
end
