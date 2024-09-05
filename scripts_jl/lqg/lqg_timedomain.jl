using multiwfs
using multiwfs: rms
using Plots
using Distributions: MvNormal, Normal
using LinearAlgebra: diag, diagm, I
using Base.GC: gc
using DSP: freq, power
using StatsBase: mean

begin
    f_loop = 1000.0
    freq_low = 0.0079
    damp_low = 1.0
    freq_high = 2.7028
    damp_high = 0.2533
    log_lf_cost = 1.7217
    log_lf_process_noise = -8.0
    log_hf_cost = -0.6807
    log_hf_process_noise = -8.0
    Av1 = A_vib(freq_high/f_loop, damp_high)
    Av2 = A_vib(freq_low/f_loop, damp_low)
    A_ar1 = [0.995 0; 1 0]
    L = A_DM(2)
    Ã = block_diag(L, A_ar1, Av1, Av2)
    C̃ = [0 -1 0 1 0 1 0 1]
    D̃ = [1 0 0 0 0 0 0 0]' 
    B = [0; 0; 1; 0; exp10(log_hf_process_noise); 0; exp10(log_lf_process_noise); 0]
    Pw = hcat(1...)
    W = B * Pw * B'
    V = hcat(1...)
    K̃ = kalman_gain(Ã, C̃, W, V)
    Vv = [0 -1 0 1 0 exp10(log_hf_cost) 0 exp10(log_lf_cost)]
    Q = Vv' * Vv
    R = zeros(1,1)
    L = lqr_gain(Ã, D̃, Q, R)
end

begin
    N = size(Ã, 1)
    W = diagm(diag(B * B'))
    W_indices = findall(diag(W) .!= 0)
    # process_noise = MvNormal(zeros(length(W_indices)), W[W_indices,W_indices])
    x, x̂, x_ol = zeros(N), zeros(N), zeros(N)
    x[W_indices] .= rand(process_noise)
    y = C̃ * x
    ikca = (I - K̃ * C̃) * Ã
    ikcd = (I - K̃ * C̃) * D̃
    nsteps = 1_000_000
    y, y_ol = zeros(nsteps), zeros(nsteps)
    for i in 1:nsteps
        x_ol = Ã * x_ol
        x = Ã * x + D̃ * L * x̂
        w = B * rand(Normal())
        x += w
        x_ol += w
        x̂ = ikca * x̂ + K̃ * C̃ * x + ikcd * L * x̂ # Kalman update
        y[i] = (C̃ * x)[1]
        y_ol[i] = (C̃ * x_ol)[1]
    end
    gc()
    plot(y_ol, label="OL, RMSE = $(round(rms(y_ol .- mean(y_ol)), digits=3))")
    plot!(y, label="CL, RMSE = $(round(rms(y .- mean(y)), digits=3))")
end

plot_psd(psd(y_ol, 1000), normalize=false)
plot_psd!(psd(y, 1000), normalize=false)

fr = freq(psd(y, 1000))
etf2 = power(psd(y, 1000)) ./ power(psd(y_ol, 1000))
etf_analytic = 1 ./ (1 .+ lqg_controller_tf(Ã, D̃, C̃, K̃, L, 1 ./ exp.(2π .* im .* fr ./ f_loop)))
plot(fr[2:end], etf2[2:end], xscale=:log10, yscale=:log10)
plot!(fr[2:end], abs2.(etf_analytic[2:end]))
