function block_diag(matrices...)
    matrices = collect(matrices)
    total_dim = 0
    for m in matrices
        @assert ndims(m) == 2
        @assert size(m, 1) == size(m, 2)
        total_dim += size(m, 1)
    end
    A = zeros(total_dim, total_dim)
    i = 1
    for m in matrices
        k = size(m, 1)
        A[i:(i+k-1), i:(i+k-1)] .= m
        i += k
    end
    A
end

export block_diag