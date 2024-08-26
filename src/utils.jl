function block_diag(matrices...)
    matrices = collect(matrices)
    total_dim = 0
    for m in matrices
        @assert ndims(m) == 2
        @assert size(m, 1) == size(m, 2)
        total_dim += size(m, 1)
    end
    types = [typeof(m).parameters[1] for m in matrices]
    if ComplexF64 in types
        t = ComplexF64
    else
        t = Float64
    end
    A = zeros(t, total_dim, total_dim)
    i = 1
    for m in matrices
        k = size(m, 1)
        A[i:(i+k-1), i:(i+k-1)] .= m
        i += k
    end
    A
end

export block_diag