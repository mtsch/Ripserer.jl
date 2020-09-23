function prog_print(progress, args...)
    progress && printstyled(stderr, args...; color=:green)
end
function prog_println(progress, args...)
    prog_print(progress, args..., '\n')
end

simplex_name(::Type{<:Simplex{1}}) = "edges"
simplex_name(::Type{<:Simplex{2}}) = "triangles"
simplex_name(::Type{<:Simplex{3}}) = "tetrahedra"
simplex_name(::Type{<:AbstractSimplex{D}}) where D = "$D-simplices"
