using Ripserer
using Ripserer: ChainElement, chain_element_type, coefficient, simplex

for F in (Mod{2}, Mod{11}, Rational{Int})
    CE = chain_element_type(Simplex, F)

    sx = Simplex{2}(3, 1.0)
    α = F(5)
    a = CE(sx, F(13))
    b = CE(sx, F(17))
    c = CE(sx, F(21))

    @test coefficient(a) == F(13)
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
