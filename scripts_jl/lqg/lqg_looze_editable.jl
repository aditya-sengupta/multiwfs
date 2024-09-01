using LinearAlgebra: I
using multiwfs

begin
    A = [0.995 0; 0 0.995]
    C = [1 1]
    D = [0 0 -1]
    n_input_history = 3
    L = zeros(n_input_history, n_input_history)
    for i in 1:(n_input_history-1)
        L[i+1,i] = 1
    end
    Ã = block_diag(L, A)
    C̃ = hcat(D, C)
    D̃ = [1 0 0 0 0]'
    W = zeros(size(Ã))
    W[n_input_history+1:end,n_input_history+1:end] = I(2)
    V = hcat(8...)
    K̃ = kalman_gain(Ã, C̃, W, V)
    Vv = [0 0 -1 1 1]
    Q = Vv' * Vv
    R = zeros(1,1)
    L = lqr_gain(Ã, D̃, Q, R)
end
