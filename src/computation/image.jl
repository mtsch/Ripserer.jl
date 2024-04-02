"""
    _max_index(n::Int, ::Val{D})

Calculate the maximum possible index of a `Simplex{D}` constructed on `n` vertices.
"""
function _max_index(subset::Int, ::Val{D}) where {D}
    acc = _binomial(subset, Val(D+1))
    for j in 2:D
        acc += _binomial(j - 2, Val(D+1))
    end
    return acc
end

"""

"""
struct ImageBoundaryMatrix{
    Implicit,
    T<:Number,
    IndexT<:Signed,
    F,
    ColIndex<:Simplex{<:Any,<:Any,IndexT},
    R<:ReducedMatrix,
    B<:Chain{T},
    C<:Chain{T},
}
    filtration::F
    reduced::R
    buffer::B
    chain::C
    columns_to_reduce::Vector{ColIndex}
    subset::IndexT
    birth_simplices::Dict{IndexT,PersistenceInterval}
    death_simplices::Dict{IndexT,PersistenceInterval}
end

function ImageBoundaryMatrix{I}(
    ::Type{T}, filtration, columns_to_reduce, subset,
    birth_simplices::Dict{IndexT}, death_simplices::Dict{IndexT}
) where {I,T,IndexT}
    if eltype(columns_to_reduce) === Any
        ColIndex = typeof(first(columns_to_reduce))
    else
        ColIndex = eltype(columns_to_reduce)
    end
    @assert index_type(ColIndex) ≡ IndexT
    BirthT = birth_type(ColIndex)
    RowIndex = Simplex{dim(ColIndex) - 1,Tuple{Bool,BirthT},IndexT}

    columns = ColIndex[]
    foreach(columns_to_reduce) do c
        push!(columns, abs(c))
    end

    if !I
        reduced = ReducedMatrix{RowIndex,T,RowIndex}()
        buffer = Chain{T,RowIndex}()
    else
        error("not implemented yet")
    end
    sizehint!(reduced, length(columns))
    chain = Chain{T,RowIndex}()

    return ImageBoundaryMatrix{
        I,T,IndexT,typeof(filtration),ColIndex,typeof(reduced),typeof(buffer),typeof(chain)
    }(
        filtration, reduced, buffer, chain, columns, IndexT(subset),
        birth_simplices, death_simplices,
    )
end

field_type(::ImageBoundaryMatrix{<:Any,T}) where {T} = T
dim(::ImageBoundaryMatrix{<:Any,<:Any,<:Any,<:Any,S}) where {S} = dim(S) - 1
ordering(::ImageBoundaryMatrix) = Base.Order.Reverse

is_implicit(::ImageBoundaryMatrix{I}) where {I} = I
is_cohomology(::ImageBoundaryMatrix) = false

# The naming here is not ideal...
function coboundary(matrix::ImageBoundaryMatrix, simplex::AbstractCell)
    return ImageBoundary(boundary(matrix.filtration, simplex), matrix.subset)
end

<<<<<<< Updated upstream
=======
append_infinite_intervals!(_, ::ImageBoundaryMatrix) = nothing

>>>>>>> Stashed changes
"""
    ImageBoundary

Wraps over `Boundary` and adds a flag to the births indicating whether the simplex belongs
to the subset or not.
"""
struct ImageBoundary{D,I,B<:Boundary{D,I}}
    boundary::B
    subset::Int
end

function Base.iterate(ib::ImageBoundary{D}, k=1) where {D}
    next = iterate(ib.boundary, k)
    if !isnothing(next)
        sx, k = next
        idx = index(sx)
        modified_sx = Simplex{D-2}(idx, (idx > _max_index(ib.subset, Val(D-2)), birth(sx)))
        return modified_sx, k
    else
        return nothing
    end
end

"""
    interval(matrix, column, pivot, cutoff, reps)

Construct a persistence interval.
"""
function interval(matrix::ImageBoundaryMatrix, column, pivot, cutoff, reps)
    if isnothing(pivot)
        # In homology, birth simplex is nothing when column is fully reduced.
        return nothing
    else
        birth_simplex, death_simplex = simplex(pivot), column
    end

    birth_time = Float64(birth(birth_simplex)[2])
    death_time = isnothing(death_simplex) ? Inf : Float64(birth(death_simplex))
<<<<<<< Updated upstream
    if death_time - birth_time > cutoff
        # Image interval only if birth simplex is also a birth simplex in the subfiltration
        # Record the intervals they are connected to.
        sub_connection = get(matrix.birth_simplices, index(birth_simplex), nothing)
        full_connection = get(matrix.death_simplices, index(death_simplex), nothing)
        if isnothing(full_connection)
            @warn "lol"
        end
        if isnothing(sub_connection) #|| isnothing(full_connection)
            return nothing
        end

        if reps
            rep = (; representative=collect_cocycle!(matrix, column, pivot))
        else
            rep = NamedTuple()
        end
        meta = (
            ;
            birth_simplex=birth_simplex,
            death_simplex=death_simplex,
            mixup_birth=sub_connection,
            mixup_death=full_connection,
            rep...,
        )
        return PersistenceInterval(birth_time, death_time, meta)
    else
        return nothing
    end
=======

    # Image interval only if birth simplex is also a birth simplex in the subfiltration
    # Record the intervals they are connected to.
    subfiltration_interval = get(matrix.birth_simplices, index(birth_simplex), nothing)
    if isnothing(subfiltration_interval)
        return nothing
    end

    if reps
        rep = (; representative=collect_cocycle!(matrix, column, pivot))
    else
        rep = NamedTuple()
    end
    meta = (
        ;
        birth_simplex=birth_simplex,
        death_simplex=death_simplex,
        mixup_birth=subfiltration_interval,
        rep...,
    )
    return PersistenceInterval(birth_time, death_time, meta)
>>>>>>> Stashed changes
end


###
### Actual computation
###
struct ImagePersistence
    subset::Int
end

<<<<<<< Updated upstream
function _ripserer(
    ip::ImagePersistence, filtration, cutoff, verbose, field, dim_max, reps, implicit
)
=======
function _sx_vertices(filtration)
    return collect(simplex(filtration, Val(0), (i,)) for i in vertices(filtration))
end

function _ripserer(
    ip::ImagePersistence, filtration, cutoff, verbose, field, dim_max, reps, implicit
)
    D = Dict{Int,PersistenceInterval}

>>>>>>> Stashed changes
    subfiltration = subset(filtration, ip.subset)

    full_result = PersistenceDiagram[]
    sub_result = PersistenceDiagram[]
    image_result = PersistenceDiagram[]

    # Subfiltration
    sub_edges = edges(subfiltration)
<<<<<<< Updated upstream
    sub_matrix = BoundaryMatrix{false}(field, subfiltration, sub_edges)
    sub_zeroth = compute_intervals!(sub_matrix, cutoff, verbose, _reps(reps, dim))
    push!(sub_result, sub_zeroth)
    birth_simplices = Dict(index(birth_simplex(i)) => i for i in sub_zeroth)

    # Full filtration
    full_edges = edges(filtration)
    full_matrix = BoundaryMatrix{false}(field, filtration, full_edges)
    full_zeroth = compute_intervals!(full_matrix, cutoff, verbose, _reps(reps, dim))
    push!(full_result, full_zeroth)
    death_simplices = Dict(
=======
    sub_matrix = BoundaryMatrix{false}(
        field, subfiltration, _sx_vertices(subfiltration), sub_edges
    )
    sub_zeroth = compute_intervals!(sub_matrix, cutoff, verbose, _reps(reps, dim))
    push!(sub_result, sub_zeroth)
    birth_simplices = D(index(birth_simplex(i)) => i for i in sub_zeroth)

    # Full filtration
    full_edges = edges(filtration)
    full_matrix = BoundaryMatrix{false}(
        field, filtration, _sx_vertices(filtration), full_edges
    )
    full_zeroth = compute_intervals!(full_matrix, cutoff, verbose, _reps(reps, dim))
    push!(full_result, full_zeroth)
    death_simplices = D(
>>>>>>> Stashed changes
        index(death_simplex(i)) => i for i in full_zeroth if !isnothing(death_simplex(i))
    )

    # Image
    image_matrix = ImageBoundaryMatrix{false}(
        field, filtration, full_edges, ip.subset,
        birth_simplices, death_simplices,
    )
    image_zeroth = compute_intervals!(
        image_matrix, cutoff, verbose, _reps(reps, dim)
    )
    push!(image_result, image_zeroth)

    if dim_max > 0
<<<<<<< Updated upstream
        sub_simplices = columns_to_reduce(subfiltration, sub_edges)
        full_simplices = columns_to_reduce(filtration, full_edges)
        for dim in 1:dim_max
            # Sub homology
            @debug "sub"
            sub_matrix = BoundaryMatrix{false}(field, subfiltration, sub_simplices)
=======
        sub_birth_candidates = sub_edges
        sub_columns = columns_to_reduce(subfiltration, sub_edges)
        full_birth_candidates = full_edges
        full_columns = columns_to_reduce(filtration, full_edges)
        for dim in 1:dim_max
            # Sub homology
            sub_matrix = BoundaryMatrix{false}(
                field, subfiltration, sub_birth_candidates, sub_columns
            )
>>>>>>> Stashed changes
            sub_intervals = compute_intervals!(
                sub_matrix, cutoff, verbose, _reps(reps, dim)
            )
            push!(sub_result, sub_intervals)
<<<<<<< Updated upstream
            birth_simplices = Dict(index(birth_simplex(i)) => i for i in sub_intervals)

            # Full homology
            @debug "full"
            full_matrix = BoundaryMatrix{false}(field, filtration, full_simplices)
=======
            birth_simplices = D(index(birth_simplex(i)) => i for i in sub_intervals)

            # Full homology
            full_matrix = BoundaryMatrix{false}(
                field, filtration, full_birth_candidates, full_columns
            )
>>>>>>> Stashed changes
            full_intervals = compute_intervals!(
                full_matrix, cutoff, verbose, _reps(reps, dim)
            )
            push!(full_result, full_intervals)
<<<<<<< Updated upstream
            death_simplices = Dict(index(death_simplex(i)) => i for i in full_intervals)

            # Image homology
            @debug "image"
            image_matrix = ImageBoundaryMatrix{false}(
                field, filtration, full_simplices, ip.subset,
=======
            death_simplices = D()#index(death_simplex(i)) => i for i in full_intervals)
            dss = [death_simplex(i) for i in full_intervals]

            # Image homology
            image_matrix = ImageBoundaryMatrix{false}(
                field, filtration, full_columns, ip.subset,
>>>>>>> Stashed changes
                birth_simplices, death_simplices,
            )
            image_intervals = compute_intervals!(
                image_matrix, cutoff, verbose, _reps(reps, dim)
            )
            push!(image_result, image_intervals)

            if dim < dim_max
<<<<<<< Updated upstream
                sub_simplices = columns_to_reduce(filtration, sub_simplices)
                full_simplices = columns_to_reduce(filtration, full_simplices)
=======
                sub_birth_candidates = sub_matrix.columns_to_reduce
                sub_columns = columns_to_reduce(filtration, sub_columns)
                full_birth_candidates = full_matrix.columns_to_reduce
                full_columns = columns_to_reduce(filtration, full_columns)
>>>>>>> Stashed changes
            end
        end
    end
    return sub_result, full_result, image_result
end
<<<<<<< Updated upstream
=======

struct InvolutedImage
    subset::Int
end

function _ripserer(
    ip::InvolutedImage, filtration, cutoff, verbose, field, dim_max, reps, implicit
)
    I = index_type(simplex_type(filtration, 1))
    death_simplices = Dict{I,PersistenceInterval}() # not needed

    subfiltration = subset(filtration, ip.subset)

    result = PersistenceDiagram[]

    # Subfiltration
    sub_zeroth, sub_to_reduce, sub_to_skip = zeroth_intervals(
        subfiltration, cutoff, verbose, field, _reps(reps, 0)
    )
    birth_simplices = Dict{I,PersistenceInterval}(
        index(birth_simplex(i)) => i for i in sub_zeroth
    )

    # Full filtration
    full_zeroth, full_to_reduce, full_to_skip = zeroth_intervals(
        filtration, -1, verbose, field, _reps(reps, 0)
    )

    # Image - need to do subfiltration as homology
    sub_matrix = BoundaryMatrix{false}(
        field, subfiltration, _sx_vertices(subfiltration), edges(subfiltration)
    )
    sub_intervals = compute_intervals!(sub_matrix, cutoff, verbose, false)
    #birth_simplices = Dict{I,PersistenceInterval}(
    #    index(birth_simplex(i)) => i for i in sub_intervals
    #)

    image_matrix = ImageBoundaryMatrix{false}(
        field, filtration, edges(filtration), ip.subset,
        birth_simplices, death_simplices,
    )
    image_zeroth = compute_intervals!(
        image_matrix, cutoff, verbose, _reps(reps, dim)
    )
    push!(result, image_zeroth)

    if dim_max > 0
        sub_comatrix = CoboundaryMatrix{true}(
            field, subfiltration, sub_to_reduce, sub_to_skip
        )
        full_comatrix = CoboundaryMatrix{true}(
            field, filtration, full_to_reduce, full_to_skip
        )
        for dim in 1:dim_max
            # Subfiltration
            sub_intervals = compute_intervals!(
                sub_comatrix, cutoff, verbose, _reps(reps, dim)
            )
            birth_simplices = Dict{I,PersistenceInterval}(
                index(birth_simplex(i)) => i for i in sub_intervals
            )

            # Get nonzero columns in full filtration
            columns, inf_births = compute_death_simplices!(full_comatrix, verbose, cutoff)

            # Image homology only using these columns.
            image_matrix = ImageBoundaryMatrix{false}(
                field, filtration, columns, ip.subset, birth_simplices, death_simplices
            )
            image_intervals = compute_intervals!(
                image_matrix, cutoff, verbose, _reps(reps, dim)
            )
            push!(result, image_intervals)

            if dim < dim_max
                sub_comatrix = next_matrix(sub_comatrix, verbose)
                full_comatrix = next_matrix(full_comatrix, verbose)
            end
        end
    end
    return result
end
>>>>>>> Stashed changes
