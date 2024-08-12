using DSP

mutable struct data
    a::Vector{Float64}
end

N = 50000
n = N ÷ 8
noverlap = 2048

d = data(rand(N))

welch_pgram(d.a, n, noverlap)
d.a = rand(N)
welch_pgram(d.a, n, noverlap)