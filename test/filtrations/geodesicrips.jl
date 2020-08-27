using LightGraphs
using Ripserer
using Test

using Ripserer: CircumSimplex, circum, dist

@testset "Circumsimplex" begin
    @testset "Constructors" begin
        @test CircumSimplex{1}([2, 1], 1) == CircumSimplex{1}(1, 1, 2)
        @test CircumSimplex{1}((1, 3), 0f0) == CircumSimplex{1}(2, 0f0, 1f0)
        @test CircumSimplex{2}([5, 2, 1], 1.0, 2.0) == CircumSimplex{2}(5, 1.0, 2.0)
    end

    @testset "AbstractSimplex interface" begin
        for D in (0, 2),
            T in (Float32, Int),
            I in (Int32, Int64)
            @testset "CircumSimplex{$D, $T, $I}" begin
                b = T(9)
                c = T(18)
                i = I(100)
                sx = CircumSimplex{D}(i, b, c)

                @testset "Type stuff" begin
                    @test sx ≡ CircumSimplex{D, T, I}(i, b, c)
                    @test typeof(sx) ≡ CircumSimplex{D, T, I}
                    D > 0 && @test_throws DomainError CircumSimplex{-D}(i, b, c)
                    @test eltype(sx) == I
                end

                @testset "Getters" begin
                    @test index(sx) == i
                    @test index(-sx) == i
                    @test birth(sx) ≡ b
                    @test circum(sx) ≡ c
                end

                @testset "Equality, hashing" begin
                    @test sx == sx
                    @test sx != CircumSimplex{D + 1}(i, b, c)
                    @test sx == CircumSimplex{D}(i, b + 1, c)
                    @test isequal(sx, CircumSimplex{D}(i, b + 1, c))
                    @test hash(sx) == hash(CircumSimplex{D}(i, b, c + 1))
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
                    @test sx < CircumSimplex{D}(I(i + 1), b + 1, c)
                    @test sx > CircumSimplex{D}(I(i - 1), b - 1, c)

                    @test sx < CircumSimplex{D}(I(i - 1), b, c)
                    @test sx > CircumSimplex{D}(I(i + 1), b, c)

                    @test sx < CircumSimplex{D}(i, b, c + 1)
                    @test sx > CircumSimplex{D}(i, b, c - 1)
                end

                @testset "Array interface, vertices" begin
                    verts = vertices(sx)

                    @test eltype(sx) == eltype(verts)
                    @test length(sx) == length(verts) == D + 1
                    @test size(sx) == (D + 1,)
                    @test firstindex(sx) == 1
                    @test lastindex(sx) == D + 1

                    @test CircumSimplex{D}(verts, birth(sx), circum(sx)) ≡ sx

                    for (i, v) in enumerate(sx)
                        @test v == verts[i]
                    end

                    @test begin @inferred vertices(sx); true end
                end

                @testset "Printing" begin
                    @test sprint(show, sx) ==
                        "+CircumSimplex{$D}($(vertices(sx)), $(b))"
                    @test sprint(show, -sx) ==
                        "-CircumSimplex{$D}($(vertices(sx)), $(b))"
                    @test sprint((i, s) -> show(i, MIME"text/plain"(), s), sx) ==
                        "$D-dimensional CircumSimplex(index=$i, birth=$b):\n  +$(vertices(sx))"
                end
            end
        end
    end
end

function torus_graph(n, m)
    adj = zeros(Int, n*m, n*m)
    cartesians = CartesianIndices((n, m))
    linears = LinearIndices((n, m))
    for vertex_cart in cartesians
        for δi in -1:1, δj in -1:1
            δi == 0 && δj == 0 && continue
            neighbor_cart = vertex_cart + CartesianIndex(δi, δj)
            neighbor_lin = linears[mod1(neighbor_cart[1], n), mod1(neighbor_cart[2], m)]
            vertex_lin = linears[vertex_cart]
            adj[vertex_lin, neighbor_lin] = 1
        end
    end
    return SimpleGraph(adj)
end

function projective_graph(n, m)
    adj = zeros(Int, n*m, n*m)
    cartesians = CartesianIndices((n, m))
    linears = LinearIndices((n, m))
    for vertex_cart in cartesians
        for δi in -1:1, δj in -1:1
            δi == 0 && δj == 0 && continue
            ni = (vertex_cart + CartesianIndex(δi, δj))[1]
            nj = (vertex_cart + CartesianIndex(δi, δj))[2]
            if 1 ≤ ni ≤ n
                ni = mod1(ni, n)
                nj = mod1(nj, m)
            else
                ni = mod1(ni, n)
                nj = mod1(-nj, m)
            end
            neighbor_lin = linears[ni, nj]
            vertex_lin = linears[vertex_cart]
            adj[vertex_lin, neighbor_lin] = 1
        end
    end
    return SimpleGraph(adj)
end

@testset "GeodesicRips" begin
    @testset "4-vertex cycle - make sure small holes don't get skipped" begin
        cycle = cycle_graph(4)
        @test ripserer(GeodesicRips(cycle))[2] == [(1, 2)]
    end

    @testset "Unlike Rips, finds equilateral triangle in 6-vertex cycle" begin
        cycle = cycle_graph(6)
        grips = GeodesicRips(cycle)
        vs = vertices(ripserer(grips)[2][1].death_simplex)
        @test vs[1] - vs[2] == 2
        @test vs[2] - vs[3] == 2
        @test vs[3] - vs[1] + 6  == 2
    end

    @testset "Bug" begin
        # Example was generated randomly. Full precision is needed to make it fail.
        cycle = cycle_graph(6)
        weights = [
            0.0 0.7740469607334781 0.0 0.0 0.0 0.23851705360578257
            0.7740469607334781 0.0 0.8975565787378181 0.0 0.0 0.0
            0.0 0.8975565787378181 0.0 0.796543749364292 0.0 0.0
            0.0 0.0 0.796543749364292 0.0 0.9905179751529374 0.0
            0.0 0.0 0.0 0.9905179751529374 0.0 0.05925349646304379
            0.23851705360578257 0.0 0.0 0.0 0.05925349646304379 0.0
        ]
        grips = GeodesicRips(cycle, weights)
        @test_broken ripserer(grips) == ripserer(dist(grips))
    end

    for m in (2, 17)
        @testset "18 cycle graph modulus=$m" begin
            grips = GeodesicRips(cycle_graph(18))
            @test ripserer(grips; dim_max=2, modulus=2) ==
                ripserer(dist(grips); dim_max=2, modulus=m)
        end
    end

    for m in (2, 5)
        @testset "12×9 torus graph modulus=$m" begin
            torus = torus_graph(12, 9)
            result = ripserer(GeodesicRips(torus), modulus=m)

            @test result == ripserer(dist(GeodesicRips(torus)), modulus=m)
            @test result[1] == vcat(fill((0, 1), 107), [(0, Inf)])
            @test result[2][1] == (1, 3)
            @test result[2][2] == (1, 4)
            @test circum(result[2][1].death_simplex) == 9
            @test circum(result[2][2].death_simplex) == 12
        end
    end

    @testset "5×7 projective plane graph modulus=2" begin
        plane = projective_graph(5, 7)
        result = ripserer(GeodesicRips(plane), dim_max=3)

        @test result == ripserer(dist(GeodesicRips(plane)), dim_max=3)
        @test length(result[2]) == 2
        @test length(result[3]) == 1
        @test circum(result[2][1].death_simplex) == 5
        @test circum(result[2][2].death_simplex) == 5
    end

    @testset "5×7 projective plane graph modulus=3" begin
        plane = projective_graph(5, 7)
        result = ripserer(GeodesicRips(plane), modulus=3, dim_max=3)

        @test result == ripserer(dist(GeodesicRips(plane)), modulus=3, dim_max=3)
        @test length(result[2]) == 1
        @test circum(result[2][1].death_simplex) == 5
    end
end
