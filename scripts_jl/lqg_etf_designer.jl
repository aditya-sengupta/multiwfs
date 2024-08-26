using ControlSystems: are, Discrete
using multiwfs
using LinearAlgebra: I
using Plots
using Distributions

begin
    f_loop = 1000.0
    Av1 = real.(A_vib(20.0/f_loop, 0.1))
    Av2 = real.(A_vib(0.001/f_loop, 0.1))
    A = block_diag(Av1, Av2)
    B = reshape([1 0 1 0], (4,1))
    C = reshape([1 0 1 0], (1,4))
    Ccost = reshape([1 0 1 0], (1,4))
    fr = exp10.(-4:0.01:log10(f_loop/2))
    s = 2π * im * fr ./ f_loop
    W, V = 1e-10 * I(4), zeros(1,1)
    Q = Ccost' * Ccost
    R = zeros(1,1)
    Pobs = real.(are(Discrete(1/f_loop), A', C', W, V))
    K = Pobs * C' * inv(C * Pobs * C' + V)
    Pcon = real.(are(Discrete(1/f_loop), A, B, Q, R))
    L = -inv(R + B' * Pcon * B) * (B' * Pcon * A)
    dstf = dynamic_system_tf.(s, Ref(A), Ref(B), Ref(C))
    lqgf = lqg_tf.(s, Ref(A), Ref(B), Ref(C), Ref(K), Ref(L))
    plot(fr, abs2.(lqgf ./ dstf), color=3, xscale=:log10, yscale=:log10, legend=nothing, xlabel="Frequency (Hz)", ylabel="|RTF|²")
end

# Next steps:
# 2. Add in the Kalman filter part so it's just an end-to-end controller TF
# 3. Put the resulting analytic ETF on sliders for the four parameters