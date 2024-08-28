using DSP: freq, power
using LinearAlgebra: I
using Plots
using multiwfs
using Distributions

σ = 1e-0
noise_dist = Normal(0, σ)
Nstep = 100_000
f_loop = 1000.0
times = 0.0:(1/f_loop):(Nstep/f_loop)
Av1 = real.(A_vib(50.0/f_loop, 0.05))
A = Av1
B = reshape([1 0], (2,1))
C = reshape([1 0], (1,2))
fr = freq(psd(rand(Nstep), f_loop))[2:end]
s = 2π * im * fr ./ f_loop
z = exp.(s)
dstf = 1e-2 ./ abs2.(1 .- A[1,1] * exp.(-s) - A[1,2] * exp.(-2s)) 
x = zeros(2)
ys = Float64[]
for _ in 1:Nstep
    global x = A*x
    global x[1] += rand(noise_dist)
    push!(ys, (C*x)[1])
end
olp = psd(ys, f_loop)
p = plot_psd(fr, power(olp)[2:end] / (σ/2)^2, normalize=false, label="OL PSD, time-domain", legend=:left, color=1)
plot!(fr, abs.(dstf), xscale=:log10, yscale=:log10, label="Analytic OL PSD", color=:black, ls=:dash)
p

