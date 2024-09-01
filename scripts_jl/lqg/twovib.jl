using Plots
using LinearAlgebra: I
using multiwfs
using Base.GC: gc

begin
    Av1 = A_vib(1.0/f_loop, 0.99)
    Av2 = A_vib(20.0/f_loop, 0.99)
    L = A_DM(2)
    Ã = block_diag(L, Av1, Av2)
    C̃ = [0 -1 0 1 0 1]
    D̃ = [1 0 0 0 0 0]' 
    B = [0; 0; 1; 0; 1e6; 0]
    Pw = hcat(1...)
    W = B * Pw * B'
    V = hcat(1...)
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
