"""
    isprime(n)

Return `true` if `n` is a prime number.
"""
function isprime(n)
    if iseven(n) || n < 2
        n == 2
    else
        p = 3
        q = n / p
        while p ≤ q
            iszero(n % p) && return false
            p += 2
            q = n / p
        end
        true
    end
end

"""
    Binomial(n_max, k_max)

Table of precomputed binomial coefficients up to `n_max` and `k_max`. Can be called like a
function and should be identical to [`binomial`](@ref) for values of `0 ≤ n ≤ n_max` and
`0 ≤ k ≤ k_max`
"""
struct Binomial
    table::Matrix{Int64}
end

function Binomial(n, k)
    table = zeros(Int, n+1, k+1)
    for i in 1:n+1
	table[i, 1] = 1;
	for j in 2:min(i, k+1)
	    table[i, j] = table[i-1, j-1] + table[i-1, j];
	    if (i <= k)
                table[i, i] = 1
            end
        end
    end
    Binomial(table)
end

Base.show(io::IO, bin::Binomial) =
    print(io, "Binomial$(size(bin.table) .- 1)")
(bin::Binomial)(n, k) =
    bin.table[n+1, k+1]

# TODO
@generated function Base.inv(value::Int, ::Val{M}) where M
    if M > 2
        isprime(M) || throw(DomainError(M, "modulus not prime"))
        inverse_arr = fill(0, M-1)
        inverse_arr[1] = 1
        for i in 2:M-1
            inverse_arr[i] = M - (inverse_arr[M % i] * floor(Int, M / i)) % M;
        end
        inverse = (inverse_arr...,)

        :($inverse[value])
    else
        :(value)
    end
end
