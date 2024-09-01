using Plots
using LinearAlgebra: I
using multiwfs
using Base.GC: gc

begin
    A_ar1 = [0.999 0; 1 0]
    Av1 = A_vib(0.001/f_loop, 0.5)
    L = A_DM(2)
    Ã = block_diag(L, A_ar1, Av1)
    C̃ = [0 -1 0 1 0 1]
    D̃ = [1 0 0 0 0 0]' 
    B = [0; 0; 1; 0; 1; 0]
    Pw = hcat(1...)
    W = B * Pw * B'
    V = hcat(0.1...)
    K̃ = kalman_gain(Ã, C̃, W, V)
    Vv = [0 -1 0 1 0 1]
    Q = Vv' * Vv
    R = zeros(1,1)
    L = lqr_gain(Ã, D̃, Q, R)
    f_loop = 1000.0
    fr = exp10.(-4:0.01:log10(f_loop/2))
    s = 2π * im * fr ./ f_loop
    z = exp.(s)
    zinvs = 1 ./ z
    gc()
    plot(fr, abs2.(1 ./ (1 .+ lqg_controller_tf(Ã, D̃, C̃, K̃, L, zinvs))), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="|ETF|²", legend=nothing)
end
