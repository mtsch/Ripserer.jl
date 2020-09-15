using CSV
using SparseArrays

function load_points(filename)
    # ripser uses Float32
    table = CSV.read(filename, header=0, type=Float32, delim=' ')
    nrow, dim = size(table)
    result = Vector{NTuple{dim, Float32}}(undef, nrow)
    for i in 1:nrow
        result[i] = Tuple(table[i, :])
    end
    return result
end

function load_dist(filename)
    return Matrix(CSV.read(filename, header=0, type=Float32, delim=' '))
end

function load_sparse(filename)
    table = CSV.read(filename, header=0, type=Float32, delim=' ')
    I = table[:, 1] .+ 1
    J = table[:, 2] .+ 1
    V = table[:, 3] .+ 1
    sparse([I; J], [J; I], [V; V])
end

function load_dipha(filename)
    open(filename, "r") do f
        magic = read(f, Int)
        @assert magic == 8067171840
        type = read(f, Int)
        @assert type == 1
        len = read(f, Int)
        dim = read(f, Int)
        if dim == 2
            m = read(f, Int)
            n = read(f, Int)
            result = Array{Float64}(undef, (m, n))
            read!(f, result)

            return result
        elseif dim == 3
            m = read(f, Int)
            n = read(f, Int)
            o = read(f, Int)
            result = Array{Float64}(undef, (m, n, o))
            read!(f, result)

            return result
        else
            error("not implemented for dim=$dim")
        end
    end
end

function load_data(filename)
    ext = splitext(filename)[2]
    if ext == ".pts"
        return load_points(filename)
    elseif ext == ".dist"
        return load_dist(filename)
    elseif ext == ".spdist"
        return load_sparse(filename)
    elseif ext == ".dipha"
        return load_dipha(filename)
    else
        error("unsupported filetype \"$(ext)\"")
    end
end
