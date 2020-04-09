using Ripserer:
    edges, is_distance_matrix, apply_threshold,
    isprime, Binomial, dist,
    set_coef, inv_mod,
    coboundary

@testset "rips" begin
    @testset "distancematrix" begin
        @testset "edges" begin
            @testset "dense" begin
                dist = [0 1 2 3;
                        1 0 4 5;
                        2 4 0 4;
                        3 5 4 0]
                res = edges(dist)
                @test res[1] == (1, (2, 1))
                @test res[2] == (2, (3, 1))
                @test res[3] == (3, (4, 1))
                @test res[4] == (4, (4, 3))
                @test res[5] == (4, (3, 2))
                @test res[6] == (5, (4, 2))
            end

            @testset "sparse" begin
                dist = sparse([0 1 0 3;
                               1 0 4 5;
                               0 4 0 4;
                               3 5 4 0])
                res = edges(dist)
                @test res[1] == (1, (2, 1))
                @test res[2] == (3, (4, 1))
                @test res[3] == (4, (4, 3))
                @test res[4] == (4, (3, 2))
                @test res[5] == (5, (4, 2))
            end

            @testset "number of edges" begin
                for i in 1:10
                    dist = rand_dist_matrix(100, 0.5)
                    n_edges = nnz(dist) ÷ 2
                    @test length(edges(dist)) == n_edges
                end
            end
        end

        @testset "is_distance_matrix" begin
            @test is_distance_matrix([0 1 0;
                                      1 0 2;
                                      0 2 0])
            @test is_distance_matrix([0 1 0 0;
                                      1 0 2 0;
                                      0 2 0 0;
                                      0 0 0 0])
            @test !is_distance_matrix([1 1 0;
                                       1 0 2;
                                       0 2 0])
            @test !is_distance_matrix([0 2 0;
                                       1 0 2;
                                       0 2 0])
            @test !is_distance_matrix([0 1 0 0;
                                       1 0 2 0;
                                       0 2 0 0])
        end

        @testset "apply_threshold" begin
            for t in 0:0.2:1
                dist_d = apply_threshold(rand_dist_matrix(10), t)
                dist_s = apply_threshold(rand_dist_matrix(10, 0.5), t)

                @test maximum(dist_d) ≤ t
                @test maximum(dist_s) ≤ t
                @test is_distance_matrix(dist_d)
                @test is_distance_matrix(dist_s)
                @test issparse(dist_d)
                @test issparse(dist_s)
            end
        end
    end

    @testset "rips" begin
        @testset "isprime" begin
            @test !isprime(1)
            @test isprime(2)
            @test isprime(3)
            @test !isprime(4)
            @test isprime(5)
            @test !isprime(6)
            @test isprime(7)
            @test !isprime(8)
            @test !isprime(9)
            @test !isprime(10)
            @test isprime(11)
        end

        @testset "Binomial" begin
            bin = Binomial(10, 15)
            for n in 1:10, k in 1:15
                @test bin(n, k) == binomial(n, k)
            end
        end

        @testset "RipsFiltration" begin
            flt = RipsFiltration{3}([0 1 2 0;
                                  1 0 3 0;
                                  2 3 0 0;
                                  0 0 0 0], 1)
            @test length(flt) == 4
            @test dist(flt, 3, 3) == 0
            @test dist(flt, 1, 2) == 1
            @test dist(flt, 1, 2) == 1
        end
    end

    @testset "Simplex" begin
        for M in (2, 17, 7487)
            for i in (1, 536, typemax(Int32))
                c = rand(Int64)
                d = rand(Float64)
                @test index(Simplex{M}(d, i, c)) == i
                @test coef(Simplex{M}(d, i, c)) == mod(c, M)
                @test diam(Simplex{M}(d, i, c)) == d
            end
        end
        @test_throws DomainError Simplex{-3}(rand(), 1, 1)
        @test_throws DomainError Simplex{7497}(rand(), 1, 1)
    end

    @testset "index(::Vector), vertices" begin
        flt = RipsFiltration{2}(rand_dist_matrix(10), 5)
        @test vertices(flt, Simplex{2}(rand(), 1, 1), 2) == [3, 2, 1]
        @test vertices(flt, Simplex{2}(rand(), 2, 1), 3) == [5, 3, 2, 1]
        @test vertices(flt, Simplex{2}(rand(), 3, 1), 1) == [3, 2]
        @test vertices(flt, Simplex{2}(rand(), 4, 1), 4) == [6, 5, 4, 2, 1]
        @test vertices(flt, Simplex{2}(rand(), 5, 1), 2) == [5, 2, 1]

        for i in 1:10
            @test index(flt, vertices(flt, Simplex{2}(rand(), i, 1), 5)) == i
        end
    end

    @testset "set_coef" begin
        @test set_coef(Simplex{3}(1, 2, 1), 2) == Simplex{3}(1, 2, 2)
        @test set_coef(Simplex{2}(2, 2, 1), 3) == Simplex{2}(2, 2, 1)
    end

    @testset "arithmetic" begin
        @test Simplex{3}(1.0, 3, 2) * 2 == Simplex{3}(1.0, 3, 1)
        @test 2 * Simplex{3}(2.0, 3, 2) == Simplex{3}(2.0, 3, 1)
        @test -Simplex{5}(10, 1, 1) == Simplex{5}(10, 1, 4)

        @test inv_mod(Val(2), 1) == 1
        @test_throws DomainError inv_mod(Val(4), 1)
        @test_throws DivideError inv_mod(Val(3), 0)

        for i in 1:16
            @test Simplex{17}(i*10, 10, i) / i == Simplex{17}(i*10, 10, 1)
            @test -Simplex{17}(i*10, 8, i) == Simplex{17}(i*10, 8, 17 - i)
        end
    end

    @testset "coboundary" begin
        @testset "line cofaces" begin
            flt = RipsFiltration{2}(Float64[0 1 3 4 5;
                                         1 0 3 4 5;
                                         3 3 0 9 9;
                                         4 4 9 0 9;
                                         5 5 9 9 0], 1)
            cb = Simplex{2, Float64}[]
            for sx in coboundary(flt, Simplex{2}(flt, 1.0, (2, 1), 1), 1)
                push!(cb, sx)
            end
            @test cb == [Simplex{2}(flt, 5.0, (5, 2, 1), 1),
                         Simplex{2}(flt, 4.0, (4, 2, 1), 1),
                         Simplex{2}(flt, 3.0, (3, 2, 1), 1)]
        end

        @testset "full graph" begin
            n = 100
            dists = ones(n, n)
            for i in 1:size(dists, 1)
                dists[i, i] = 0
            end
            flt = RipsFiltration{3}(dists, 5)

            for dim in 1:5
                cob = Simplex{3, Float64}[]
                sx = Simplex{3}(1.0, 10, 1)
                for c in coboundary(flt, sx, dim)
                    push!(cob, sx)
                end
                sx_vxs = vertices(flt, sx, dim)
                all(x -> issubset(sx_vxs, vertices(flt, x, dim+1)), cob)
                @test length(cob) == n - dim - 1
            end
        end

        @testset "icosahedron" begin
            dim = 2
            flt = RipsFiltration{7}(icosahedron, 4)

            sx_vxs = (10, 3, 1)
            sx = Simplex{7}(flt, diam(flt, sx_vxs), sx_vxs, 1)
            for cf in coboundary(flt, sx, dim)
                cf_vxs = vertices(flt, cf, dim+1)
                @test diam(cf) == diam(flt, cf_vxs)
                @test issubset(sx_vxs, cf_vxs)
            end
        end

        @testset "torus" begin
            dim = 2
            flt = RipsFiltration{3}(torus(16), 4)

            sx_vxs = (16, 8, 1)
            sx = Simplex{3}(flt, diam(flt, sx_vxs), sx_vxs, 1)
            for cf in coboundary(flt, sx, dim)
                cf_vxs = vertices(flt, cf, dim+1)
                @test diam(cf) == diam(flt, cf_vxs)
                @test issubset(sx_vxs, cf_vxs)
            end
        end
    end
end
