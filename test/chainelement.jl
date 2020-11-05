using Test
using Ripserer

using Ripserer:
    ChainElement,
    PackedElement,
    chain_element_type,
    simplex,
    coefficient,
    index_overflow_check

for S in (Simplex{2,Float64,Int}, Simplex{3,Float64,Int32})
    for F in (Mod{2}, Mod{11}, Rational{Int}, Mod{251}, Mod{257})
        @testset "ChainElement with $S, $F" begin
            simplex_1 = S(5, 3.0)
            simplex_2 = S(3, 5.0)

            coef_1 = F(13)
            coef_2 = F(17)
            coef_3 = F(7919)

            @testset "Constructing a chain element is inferrable" begin
                @test begin
                    @inferred chain_element_type(S, F)
                    true
                end
                @test begin
                    @inferred chain_element_type(simplex_1, coef_1)
                    true
                end
                @test chain_element_type(S, F) === chain_element_type(simplex_1, coef_1)
            end

            Element = chain_element_type(S, F)
            @testset "Printing" begin
                @test sprint(show, Element(simplex_1)) == "$(simplex_1) => $(one(F))"
                @test sprint(show, Element(simplex_2, coef_1)) ==
                      "$(simplex_2) => $(coef_1)"
            end
            @testset "Simplex, coefficient, birth, index and iteration" begin
                a = Element(-simplex_1, coef_3)
                s, c = a
                @test s === simplex_1 === simplex(a)
                @test c === -coef_3 === coefficient(a)
                @test first(a) === s
                @test last(a) === c
                @test a[1] === s
                @test a[2] === c
                @test length(a) == 2
                @test firstindex(a) == 1
                @test lastindex(a) == 2
                @test first(a) === s
                @test last(a) === -coef_3
                @test eltype(a) == Union{S,F}
                @test eltype(typeof(a)) == Union{S,F}
                @test convert(S, a) == simplex(a)
                @test_throws BoundsError a[3]
                @test_throws BoundsError a[0]

                @test birth(a) == birth(s)
                @test index(a) == index(s)
            end
            @testset "Arithmetic" begin
                α = F(5)
                a = @inferred Element(-simplex_1, coef_1)
                b = @inferred Element(simplex_1, coef_2)
                c = @inferred Element(simplex_1, coef_3)

                @test (a + b) + c == a + (b + c)
                @test a + b == b + a
                @test α * b == b * α
                @test a - a == zero(a)
                @test iszero(a - a)
                @test b / coefficient(b) == oneunit(a)
                @test b + zero(a) == b
                @test c * one(a) == c
                @test a - b == a + -b
                @test a / α == a * inv(α)
                @test α * (b + c) == α * b + α * c

                @test one(a) == one(typeof(a)) == one(α) == one(typeof(α))

                @test Element(simplex_1) === oneunit(Element(simplex_1))
                @test Element(-simplex_1) === -Element(simplex_1)
                @test iszero(Element(-simplex_1) + Element(simplex_1))

                @test sign(a) == sign(-coef_1)
                @test sign(b) == sign(coef_2)

                @test_throws ArgumentError Element(simplex_2) + Element(simplex_1)
            end
            @testset "Hash, equality and order" begin
                a = @inferred Element(simplex_1, coef_1)
                b = @inferred Element(simplex_1, coef_2)
                c = @inferred Element(simplex_2, coef_1)
                d = @inferred Element(simplex_2, coef_2)

                @test a == b
                @test c == d
                @test a < c
                @test b != d

                @test hash(a) == hash(b)
                @test hash(c) == hash(d)
            end
        end
    end

    @testset "PackedElement" begin
        @testset "Construct a `PackedElement` for Simplex if n_bits < 8" begin
            for S in (Simplex{2,Float64,Int}, Simplex{3,Int,Int128})
                @test @inferred(chain_element_type(S, Mod{2})) <: PackedElement{S,Mod{2}}
                @test @inferred(chain_element_type(S, Mod{251})) <:
                      PackedElement{S,Mod{251}}
                @test @inferred(chain_element_type(S, Mod{257})) == ChainElement{S,Mod{257}}
                @test @inferred(chain_element_type(S, Rational{Int})) ==
                      ChainElement{S,Rational{Int}}
            end
        end
    end

    @testset "Overflow" begin
        @test begin
            index_overflow_check(Simplex{2,Float64,Int}, Mod{2}, 1000)
            true
        end
        @test begin
            index_overflow_check(Cube{3,Float32,4}, Mod{2}, 2_000_000_000)
            true
        end
        @test_throws OverflowError index_overflow_check(
            Simplex{5,Float64,Int}, Mod{2}, 10000
        )

        big_index = typemax(Int) >> 7
        n = Ripserer._vertices(big_index, Val(3))[1]
        @test begin
            index_overflow_check(Simplex{2,Int,Int}, Mod{2}, n)
            true
        end
        @test_throws OverflowError index_overflow_check(Simplex{2,Int,Int}, Mod{251}, n)
        # mod 257 doesn't pack, so there is no overflow.
        @test begin
            index_overflow_check(Simplex{2,Int,Int}, Mod{257}, n)
            true
        end
    end
end
