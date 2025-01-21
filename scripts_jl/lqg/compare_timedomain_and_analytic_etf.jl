using multiwfs
using multiwfs: rms
using Distributions
using Plots
using LinearAlgebra: I
using Base.GC: gc
using NPZ

begin
    f_loop = 1000.0
    A_ar1 = [0.995 0; 1 0]
    L = A_DM(2)
    Ã = block_diag(L, A_ar1)
    C̃ = [0 -1 0 1]
    D̃ = [1 0 0 0]' 
    B = [0; 0; 1; 0]
    Pw = hcat(1...)
    W = B * Pw * B'
    V = hcat(1...)
    K̃ = kalman_gain(Ã, C̃, W, V)
    Vv = [0 -1 0 1]
    Q = Vv' * Vv
    R = zeros(1,1)
    L = lqr_gain(Ã, D̃, Q, R)
    lqg = LQG(Ã, D̃, C̃, K̃, L)
end

npzwrite(
    "data/lqgfirst_lowercutoff.npz",
    Dict(
        "K" => K̃,
        "G" => L,
        "ikcA" => lqg.ikcA,
        "ikcD" => lqg.ikcD
    )
)

f_loop = 1000.0 # Hz
fr = exp10.(-4:0.01:log10(f_loop/2))

tf_analytic = transfer_function.(Ref(lqg), 2π .* im .* fr ./ f_loop)
etf_analytic = 1 ./ (1 .+ tf_analytic)
plot(fr, abs2.(etf_analytic), xscale=:log10, yscale=:log10, legend=nothing)

begin
    N = size(Ã, 1)
    x, x̂, x_ol = zeros(N), zeros(N), zeros(N)
    x += B * rand(Normal())
    y = C̃ * x
    ikca = (I - K̃ * C̃) * Ã
    ikcd = (I - K̃ * C̃) * D̃
    nsteps = 50_000
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

fr, y_psd = genpsd(y, 1000)
fr, yol_psd = genpsd(y_ol, 1000)
etf2 = y_psd ./ yol_psd
etf_analytic = 1 ./ (1 .+ transfer_function.(Ref(lqg), 2π .* im .* fr ./ f_loop))
plot(fr, etf2, xscale=:log10, yscale=:log10, label="Time-domain ETF", xlabel="Frequency (Hz)", ylabel="|ETF|²", legend=:bottomright)
plot!(fr, abs2.(etf_analytic), label="Analytic ETF")

# yay

# now let's confirm this whole thing still holds with my LQG-first HPF controller

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
    lqg = LQG(Ã, D̃, C̃, K̃, L)
end

npzwrite("data/lqgfirst_endsummer2024.npz",Dict(
    "A" => lqg.A,
    "D" => lqg.D,
    "C" => lqg.C,
    "K" => K̃,
    "G" => L,
    "ikcA" => lqg.ikcA,
    "ikcD" => lqg.ikcD
))

begin
    N = size(Ã, 1)
    x, x̂, x_ol = zeros(N), zeros(N), zeros(N)
    x += B * rand(Normal())
    y = C̃ * x
    ikca = (I - K̃ * C̃) * Ã
    ikcd = (I - K̃ * C̃) * D̃
    nsteps = 100_000
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

fr, y_psd = genpsd(y, 1000)
fr, yol_psd = genpsd(y_ol, 1000)
etf2 = y_psd ./ yol_psd
etf_analytic = 1 ./ (1 .+ transfer_function.(Ref(lqg), 2π .* im .* fr ./ f_loop))
plot(fr, etf2, xscale=:log10, yscale=:log10, label="Time-domain ETF", xlabel="Frequency (Hz)", ylabel="|ETF|²", legend=:bottomright)
plot!(fr, abs2.(etf_analytic), label="Analytic ETF", ls=:dash)

# yay yay!

