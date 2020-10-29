using Test
using DataStructures
using Ripserer

using Ripserer: Chain, ChainElement, heapmove!, clean!

for F in (Mod{2}, Mod{3}, Mod{257}, Rational{Int})
    @testset "$F" begin
        S = Simplex{1,Int,Int}
        simplices = [Simplex{1}(1,1), Simplex{1}(2,2)]
        elements = ChainElement{S,F}.(simplices)

        @testset "Constructors" begin
            chain1 = Chain{F,S}(simplices)
            chain2 = Chain{F,S}(elements)
            @test chain1 == chain2

            chain3 = Chain{F,S}()
            append!(chain3, simplices)
            chain4 = Chain{F,S}()
            append!(chain4, elements)

            @test chain1 == chain3 == chain4
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
            elements = ChainElement{S,F}.(
                S.([1, -7, 2, 3, 4, 7, 5, 6, -1], [1, 7, 1, 1, 4, 7, 5, 6, 1])
            )
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
