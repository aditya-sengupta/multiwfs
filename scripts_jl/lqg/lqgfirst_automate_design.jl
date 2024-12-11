using multiwfs
using Plots
using Optim: optimize, minimizer, Options, Newton, LBFGS

# first, we need a model with lots of free parameters
# let's say three ARs and three vibrations, totally arbitrarily

f_loop = 1000.0
fr = exp10.(-4:0.01:log10(f_loop/2))

function etf_from_optpars(pars)
    try
        a, f, d, logp, logc = pars[1:3], pars[4:6], pars[7:9], pars[10:15], pars[16:21]
        Ã = block_diag(A_DM(2), [a[1] 0; 1 0], [a[2] 0; 1 0], [a[3] 0; 1 0], A_vib(f[1]/f_loop, d[1]), A_vib(f[2]/f_loop, d[2]), A_vib(f[3]/f_loop, d[3]))
        C̃ = [0 -1 0 1 0 1 0 1 0 1 0 1 0 1]
        D̃ = [1 0 0 0 0 0 0 0 0 0 0 0 0 0]' 
        B = [0; 0; exp10(logp[1]); 0; exp10(logp[2]); 0; exp10(logp[3]); 0; exp10(logp[4]); 0; exp10(logp[5]); 0; exp10(logp[6]); 0]
        Pw = hcat(1...)
        W = B * Pw * B'
        V = hcat(1...)
        K̃ = kalman_gain(Ã, C̃, W, V)
        Vv = [0 -1 0 exp10(logc[1]) 0 exp10(logc[2]) 0 exp10(logc[3]) 0 exp10(logc[4]) 0 exp10(logc[5]) 0 exp10(logc[6])]
        Q = Vv' * Vv
        R = zeros(1,1)
        L = lqr_gain(Ã, D̃, Q, R)
        lqg = LQG(Ã, D̃, C̃, K̃, L)
        fr = exp10.(-4:0.01:log10(f_loop/2))
        tf_analytic = transfer_function.(Ref(lqg), 2π .* im .* fr ./ f_loop)
        return 1 ./ (1 .+ tf_analytic)
    catch e
        return Inf * fr
    end
end

p0 = [0.995, 0.99, 0.95, 0.1, 1.0, 10.0, 0.1, 0.5, 1.0, -8, -8, -8, 0, 0, 0, 0, 0, 0, 0, 0, 0]

etf_from_optpars(p0)

function etf_cost(etf, fr_cutoff_idx, λ=100)
    # strong regularization on the ETF = 1 at fr < fr_cutoff constraint
    lowfreq_passing_cost = λ .* sum(abs2, exp10.(abs.(abs.(etf[1:fr_cutoff_idx]) .- 1)))
    # and otherwise try to minimize the ETF
    max_rejection_cost = sum(abs2, etf[fr_cutoff_idx:end])
    return lowfreq_passing_cost + max_rejection_cost
end

# pick a cutoff frequency
fr_cutoff = 1e-3 # Hz
fr_cutoff_idx = findfirst(fr .> fr_cutoff)

# now, let's try and minimize the ETF cost?

res = optimize(pars -> etf_cost(etf_from_optpars(pars), fr_cutoff_idx), p0, Newton())

plot(fr, abs2.(etf_from_optpars(res.minimizer)), xscale=:log10, yscale=:log10)