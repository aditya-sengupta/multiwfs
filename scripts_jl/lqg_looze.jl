using Plots
using LinearAlgebra: I
using multiwfs: block_diag
using ControlSystems: are, Discrete

begin
    A = [0.995 0 0 0; 0 0.968-0.156im 0 0; 1 1 0 0; 0 0 1 0]
    B = [1 0; 0 0.5; 0 0; 0 0]
    δ = 0.654
    C = [0 0 1-δ δ]
    D = [0 δ-1 -δ]
    n_input_history = length(D)
    L = zeros(n_input_history, n_input_history)
    for i in 1:(n_input_history-1)
        L[i+1,i] = 1
    end
    Ã = block_diag(L, A)
    B̃ = vcat(zeros(3,2), B)
    C̃ = hcat(D, C)
    W = zeros(ComplexF64, size(Ã))
    W[n_input_history+1:end,n_input_history+1:end] = I(4)
    V = hcat(8...)
    K̃ = kalman_gain(Ã, C̃, W, V)
    @assert all(K̃[1:3] .+ 1 .≈ 1)
    K̃[1:3] .= 0
    Vv = [0 0 -1 0 0 1]
    Q = Vv' * Vv
    R = ones(n_input_history, n_input_history)
    L = lqr_gain(Ã, B̃, Q, R)
end
