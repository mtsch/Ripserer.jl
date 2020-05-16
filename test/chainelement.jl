using Ripserer
using Ripserer: ChainElement, PackedElement, chain_element_type, coefficient, simplex

@testset "arithmetic" begin
    for F in (Mod{2}, Mod{11}, Rational{Int}, Mod{251}, Mod{257})
        CE = @inferred chain_element_type(Simplex, F)

        sx = Simplex{2}(3, 1.0)
        α = F(5)
        a = @inferred CE(-sx, F(13))
        b = @inferred CE(sx, 17)
        c = @inferred CE(sx, 21)

        @test coefficient(a) == F(-13)
        @test coefficient(b) == F(17)
        @test coefficient(c) == F(21)
        @test simplex(a) == sx
        @test simplex(b) == sx
        @test simplex(c) == sx

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

        @test CE(sx) == oneunit(CE(sx))
        @test CE(-sx) == -CE(sx)
        @test iszero(CE(-sx) + CE(sx))
    end
end

@testset "chain_element_type" begin
    for S in (Simplex{2, Float64, Int}, Simplex{3, Int, Int128})
        @test chain_element_type(S, Mod{2}) <: PackedElement{S, Mod{2}}
        @test chain_element_type(S, Mod{2}) isa DataType
        @test chain_element_type(S, Mod{251}) <: PackedElement{S, Mod{251}}
        @test chain_element_type(S, Mod{251}) isa DataType
        @test chain_element_type(S, Mod{257}) <: ChainElement{S, Mod{257}}
        @test chain_element_type(S, Mod{257}) isa DataType
        @test chain_element_type(S, Rational{Int}) <: ChainElement{S, Rational{Int}}
        @test chain_element_type(S, Rational{Int}) isa DataType
    end
end
