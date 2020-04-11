using Ripserer:
    edges, is_distance_matrix,
    isprime, n_bits, Binomial,
    set_coef, inv_mod,
    coboundary

@testset "rips" begin
    @testset "distancematrix" begin
        @testset "edges dense" begin
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
        @testset "edges sparse" begin
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
        @testset "n edges dense" begin
            dist = rand_dist_matrix(100)
            n_edges = binomial(size(dist, 1), 2)
            @test length(edges(dist)) == n_edges
        end
        @testset "n edges sparse" begin
            for _ in 1:10
                dist = rand_dist_matrix(100, 0.5)
                n_edges = nnz(dist) ÷ 2
                @test length(edges(dist)) == n_edges
            end
        end
        @testset "is_distance_matrix" begin
            for f in (identity, sparse)
                @test is_distance_matrix(f([0 1 0;
                                            1 0 2;
                                            0 2 0]))
                @test is_distance_matrix(f([0 1 0 0;
                                            1 0 2 0;
                                            0 2 0 0;
                                            0 0 0 0]))
                @test !is_distance_matrix(f([1 1 0;
                                             1 0 2;
                                             0 2 0]))
                @test !is_distance_matrix(f([0 2 0;
                                             1 0 2;
                                             0 2 0]))
                @test !is_distance_matrix(f([0 1 0 0;
                                             1 0 2 0;
                                             0 2 0 0]))
            end
        end
    end

    @testset "helpers" begin
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
        @testset "n_bits" begin
            @test n_bits(2) == 1
            @test n_bits(3) == 2
            @test n_bits(7) == 3
            @test n_bits(Int64(typemax(UInt32))) == 32
        end
        @testset "Binomial" begin
            bin = Binomial(10, 15)
            for n in 1:10, k in 1:15
                @test bin(n, k) == binomial(n, k)
            end
        end
    end

    @testset "Simplex" begin
        @testset "index, diam, coef, set_coef" begin
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
                @test coef(Simplex{17}(1.0, 1, i) - Simplex{17}(1.0, 1, 17-i)) ==
                    coef(Simplex{17}(1.0, 1, i) + -Simplex{17}(1.0, 1, 17-i))
            end
        end
    end

    for Filtration in (RipsFiltration, SparseRipsFiltration)
        typename = string(Filtration)
        @testset "$typename" begin
            @testset "length, dist, threshold" begin
                flt = Filtration([0 1 2 9;
                                  1 0 3 9;
                                  2 3 0 4;
                                  9 9 4 0], modulus=3)
                @test length(flt) == 4
                @test dist(flt, 3, 3) == 0
                @test dist(flt, 1, 2) == 1
                @test dist(flt, 1, 3) == 2
                @test dist(flt, 3, 2) == 3
                @test threshold(flt) == 4

                flt = Filtration([0 1 2;
                                  1 0 3;
                                  2 3 0], threshold=2)
                @test length(flt) == 3
                @test dist(flt, 3, 3) == 0
                @test dist(flt, 1, 2) == 1
                @test dist(flt, 1, 3) == 2
                @test dist(flt, 3, 2) == (issparse(Filtration) ? typemax(Int) : 3)
                @test threshold(flt) == 2
            end
            @testset "index(::Vector), vertices" begin
                flt = Filtration(rand_dist_matrix(10), dim_max=5, modulus=2)
                @test vertices(flt, Simplex{2}(rand(), 1, 1), 2) == [3, 2, 1]
                @test vertices(flt, Simplex{2}(rand(), 2, 1), 3) == [5, 3, 2, 1]
                @test vertices(flt, Simplex{2}(rand(), 3, 1), 1) == [3, 2]
                @test vertices(flt, Simplex{2}(rand(), 4, 1), 4) == [6, 5, 4, 2, 1]
                @test vertices(flt, Simplex{2}(rand(), 5, 1), 2) == [5, 2, 1]

                for i in 1:10
                    @test index(flt, vertices(flt, Simplex{2}(rand(), i, 1), 5)) == i
                end
            end
            @testset "coboundary" begin
                @testset "line cofaces" begin
                    flt = Filtration(Float64[0 1 3 4 5;
                                             1 0 3 4 5;
                                             3 3 0 9 9;
                                             4 4 9 0 9;
                                             5 5 9 9 0])
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
                    flt = Filtration(dists, dim_max=5, modulus=3)

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
                    flt = Filtration(icosahedron, dim_max=4, modulus=7)

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
                    flt = Filtration(torus(16), dim_max=4, modulus=3)

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
    end
end
