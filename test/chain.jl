using Test
using DataStructures
using Ripserer

using Ripserer:
    Chain,
    ChainElement,
    PackedElement,
    chain_element_type,
    simplex,
    coefficient,
    index_overflow_check,
    heapmove!,
    clean!

@testset "ChainElements" for S in (Simplex{2,Float64,Int}, Simplex{3,Float64,Int32})
    @testset "Integers and floats are not allowed as field types" begin
        @test_throws ErrorException ChainElement{S, Int}(S(5, 3.0))
        @test_throws ErrorException ChainElement{S, UInt8}(S(5, 3.0))
        @test_throws ErrorException ChainElement{S, Float64}(S(5, 3.0))
        @test_throws ErrorException ChainElement{S, BigFloat}(S(5, 3.0))
    end

    for F in (Mod{2}, Mod{11}, Rational{Int}, Mod{251}, Mod{257})
        @testset "chain element with $S, $F" begin
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

@testset "Chain" for F in (Mod{2}, Mod{3}, Mod{257}, Rational{Int})
    @testset "$F" begin
        S = Simplex{1,Int,Int}
        simplices = [Simplex{1}(1, 1), Simplex{1}(2, 2)]
        elements = ChainElement{S,F}.(simplices)

        @testset "Constructors, basics" begin
            @test begin
                @inferred Chain{F,S}(simplices)
                @inferred Chain{F,S}(elements)
                @inferred Chain{F,S}()
                @inferred Chain{F,S}()
                true
            end
            chain1 = Chain{F,S}(simplices)
            chain2 = Chain{F,S}(elements)
            chain3 = Chain{F,S}()
            chain4 = Chain{F,S}()
            @test chain1 == chain2

            append!(chain3, simplices)
            append!(chain4, elements)

            @test chain1 == chain3 == chain4

            @test sprint(summary, chain1) == "2-element Chain{$F,$S}"

            @test simplex_type(chain1) ≡ S
            @test field_type(chain2) ≡ F
        end

        @testset "Array stuff" begin
            chain = Chain{F,S}()
            @test isempty(chain)

            push!(chain, simplices[1])
            @test length(chain) == 1
            chain[1] === elements[1]

            push!(chain, simplices[2])
            @test size(chain) == (2,)
            @test pop!(chain) === elements[2]

            empty!(chain)
            @test isempty(chain)

            resize!(chain, 2)
            chain .= simplices
            @test chain[2] === elements[2]

            @test chain[1:2] isa Chain
            @test chain[begin:end] == chain

            @test view(chain, 1:2) == chain
        end

        @testset "Heap stuff" begin
            elements =
                ChainElement{
                    S,F
                }.(S.([1, -7, 2, 3, 4, 7, 5, 6, -1], [1, 7, 1, 1, 4, 7, 5, 6, 1]))
            uniques = ChainElement{S,F}.(S.([7, 2, 3, 4, 5, 6], [7, 1, 1, 4, 5, 6]))

            fwd = Base.Order.Forward
            rev = Base.Order.Reverse

            @testset "A fresh Chain `heappop!`s nothing" begin
                chain = Chain{F,S}(fwd)
                @test heappop!(chain) ≡ nothing
            end

            @testset "pushing elements and popping finds the lowest simplex (fwd)" begin
                chain = Chain{F,S}(fwd)
                for e in elements
                    heappush!(chain, e)
                end

                @test heappop!(chain) ≡ ChainElement{S,F}(S(3, 1))
            end

            @testset "`append!`, `heapify!`, and `heappop!` (rev)" begin
                chain = Chain{F,S}(rev)
                append!(chain, elements)
                heapify!(chain)

                @test heappop!(chain) ≡ ChainElement{S,F}(S(6, 6))
            end

            @testset "adding inverses removes elements (fwd)" begin
                chain = Chain{F,S}(fwd)
                append!(chain, uniques)
                heapify!(chain)
                for e in uniques[1:2:end]
                    heappush!(chain, -e)
                end

                chain′ = copy(chain)
                for e in sort(uniques[2:2:end]; order=fwd)
                    @test heappop!(chain) == e
                end
                @test heappop!(chain) ≡ nothing
                @test heapmove!(chain′) == sort(uniques[2:2:end]; order=fwd)
                @test isempty(chain)
                @test isempty(chain′)
            end

            @testset "`clean!`" begin
                chain_fwd = Chain{F,S}(elements, fwd)
                chain_rev = Chain{F,S}(elements, rev)

                clean!(chain_fwd)
                clean!(chain_rev)

                @test chain_fwd == reverse!(chain_rev)
                @test allunique(chain_fwd)
                @test all(!iszero, chain_fwd)
            end
        end
    end
end
