using ControlSystems: are, Discrete
using DSP: freq, power
using multiwfs
using multiwfs: block_diag
using LinearAlgebra: I
using Plots

# missing something very obvious here...

begin
    f_loop = 1000.0
    Av1 = A_vib(50.0/f_loop, 0.1)
    A = block_diag(Av1, zeros(1,1))
    B = reshape([0 0 1], (3,1))
    C = reshape([1 0 -1], (1,3))
    fr = 0.16:0.16:f_loop/2
    s = 2π * im * fr
    W, V = 1e-8 * I(3), zeros(1,1)
    Q = C' * C
    R = 1e-8 * I(1)
    Pobs = are(Discrete(1/f_loop), A', C', W, V)
    K = Pobs * C' * inv(C * Pobs * C' + V)
    Pcon = are(Discrete(1/f_loop), A, B, Q, R)
    L = -inv(R + B' * Pcon * B) * (B' * Pcon * A)
    dstf = dynamic_system_tf.(s, Ref(A), Ref(B), Ref(C))
    lqgf = lqg_tf.(s, Ref(A), Ref(B), Ref(C), Ref(K), Ref(L))
end;

plot(
    fr ./ f_loop,
    f -> 1 / abs2(1 - A[1,1] * exp(-2π * im * f) - A[1,2] * exp(-2π * im * 2 * f)),
    xscale=:log10, yscale=:log10
)

begin
    x, xcon = zeros(3), zeros(3)
    x[1], xcon[1] = 1.0, 1.0
    x[2], xcon[2] = -0.0, 0.0
    ys, ycons = Float64[], Float64[]
    for _ in 1:50_000
        x = A*x
        xcon = (A+B*L)*xcon
        push!(ys, (C*x)[1])
        push!(ycons, (C*xcon)[1])
    end
    olp = psd(ys, f_loop)
    clp = psd(ycons, f_loop)
    plot_psd(freq(olp)[2:end], power(olp)[2:end], normalize=false, label="OL PSD, time-domain", legend=:bottomright, color=1)
    plot!(fr, abs2.(dstf), xscale=:log10, yscale=:log10, label="Analytic OL PSD", color=:black, ls=:dash)
    #plot_psd!(freq(clp)[2:end], power(clp)[2:end], normalize=false, label="CL PSD, time-domain", color=2)
    #plot!(fr, abs2.(lqgf) , xscale=:log10, yscale=:log10, label="Analytic CL PSD", color=2, ls=:dash)
    #plot!(fr, (power(clp) ./ power(olp))[2:end], label="RTF, time-domain", color=3)
    #plot!(fr, abs2.(lqgf ./ dstf), label="Analytic RTF", color=3, ls=:dash)
end
