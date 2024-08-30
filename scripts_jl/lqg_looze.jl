using Plots
using LinearAlgebra: I
using multiwfs: block_diag
using ControlSystems: are, Discrete

begin
    A = [0.99 0; 1 0]
    δ = 1
    B = [0; 0; 1; 0]
    C = [1-δ δ]
    D = [δ-1 -δ]
    n_input_history = length(D)
    L = zeros(n_input_history, n_input_history)
    for i in 1:(n_input_history-1)
        L[i+1,i] = 1
    end
    Ã = block_diag(L, A)
    C̃ = hcat(D, C)
    D̃ = [1 0 0 0]'
    Pw = hcat(1...)
    W = B * Pw * B'
    V = hcat(0.1...)
    K̃ = kalman_gain(Ã, C̃, W, V)
    Vv = [0 -1 0 1]
    Q = Vv' * Vv
    R = zeros(1,1)
    L = lqr_gain(Ã, D̃, Q, R)
    f_loop = 1000.0
    fr = exp10.(-4:0.01:log10(f_loop/2))
    s = 2π * im * fr ./ f_loop
    z = exp.(s)
    zinvs = 1 ./ z

    plot(fr, abs2.(1 ./ (1 .+ lqg_controller_tf(Ã, D̃, C̃, K̃, L, zinvs))), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="|ETF|²", legend=nothing, title="Delay = 2")
end

function lqg_controller_tf(A, D, C, K, G, zinvs)
    ikcA = (I - K * C) * A
    ikcD = (I - K * C) * D
    numerator = [(G * inv(I - ikcA * zinv))[1,1] for zinv in zinvs]
	denominator = [(I - G * inv(I - ikcA * zinv) * ikcD * zinv)[1,1] for zinv in zinvs]
    return numerator .* (zinvs .^ 2) ./ denominator
end
