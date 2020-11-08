using Test
using DataStructures
using Ripserer

using Ripserer: Chain, ChainElement, coefficient_type, simplex_type, heapmove!, clean!

for F in (Mod{2}, Mod{3}, Mod{257}, Rational{Int})
    @testset "Chain with $F" begin
        S = Simplex{1,Int,Int}
        simplices = [Simplex{1}(1, 1), Simplex{1}(2, 2)]
        elements = ChainElement{S,F}.(simplices)

        @testset "Constructors, basics" begin
            @test begin
                @inferred Chain{F,S}(simplices)
                @inferred Chain{F,S}(elements)
                @inferred Chain(view(elements, 1:2))
                @inferred Chain{F,S}()
                @inferred Chain{F,S}()
                true
            end
            chain1 = Chain{F,S}(simplices)
            chain2 = Chain{F,S}(elements)
            chain3 = Chain(view(elements, 1:2))
            chain4 = Chain{F,S}()
            chain5 = Chain{F,S}()

            append!(chain4, simplices)
            append!(chain5, elements)

            @test chain1 == chain2 == chain3 == chain4 == chain5

            @test sprint(summary, chain1) == "2-element Chain{$F,$S}"

            @test simplex_type(chain1) ≡ S
            @test coefficient_type(chain2) ≡ F
        end

        @testset "Internal eltypes" begin
            chain = Chain{F,S}(simplices)

            if F === Mod{2}
                @test eltype(chain.elements) <: Simplex
            elseif F === Mod{3}
                @test eltype(chain.elements) <: Ripserer.PackedElement
            else
                @test eltype(chain.elements) <: ChainElement
            end
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
            @test chain[firstindex(chain):end] == chain

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
                chain = Chain{F,S}()
                @test heappop!(chain, fwd) ≡ nothing
            end

            @testset "pushing elements and popping finds the lowest simplex (fwd)" begin
                chain = Chain{F,S}()
                for e in elements
                    heappush!(chain, e, fwd)
                end

                @test heappop!(chain, fwd) ≡ ChainElement{S,F}(S(3, 1))
            end

            @testset "`append!`, `heapify!`, and `heappop!` (rev)" begin
                chain = Chain{F,S}()
                append!(chain, elements)
                heapify!(chain, rev)

                @test heappop!(chain, rev) ≡ ChainElement{S,F}(S(6, 6))
            end

            @testset "adding inverses removes elements (fwd)" begin
                chain = Chain{F,S}()
                append!(chain, uniques)
                heapify!(chain, fwd)
                for e in uniques[1:2:end]
                    heappush!(chain, -e, fwd)
                end

                chain′ = copy(chain)
                for e in sort(uniques[2:2:end]; order=fwd)
                    @test heappop!(chain) == e
                end
                @test heappop!(chain, fwd) ≡ nothing
                @test heapmove!(chain′, fwd) == sort(uniques[2:2:end]; order=fwd)
                @test isempty(chain)
                @test isempty(chain′)
            end

            @testset "`clean!`" begin
                chain_fwd = Chain{F,S}(elements)
                chain_rev = Chain{F,S}(elements)

                clean!(chain_fwd, fwd)
                clean!(chain_rev, rev)

                @test chain_fwd == reverse!(chain_rev)
                @test allunique(chain_fwd)
                @test all(!iszero, chain_fwd)
            end
        end
    end
end
