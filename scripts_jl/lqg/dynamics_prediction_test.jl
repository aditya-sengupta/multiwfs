using multiwfs
using multiwfs: psd_ar, turbulence_ar2_loss
using StatsBase: mean
using Plots
using Distributions

function rms(x)
    x = x .- mean(x)
    return sqrt(mean(x .^ 2))
end

f_loop = 200
f_unitless = 1e-3:1e-3:0.5
fr = f_unitless .* f_loop

begin
    psd_vk = psd_von_karman(0.0)
    res = Optim.optimize(ϕ -> turbulence_ar2_loss(ϕ, psd_vk, f_unitless), zeros(2), NelderMead())
    ϕ = Optim.minimizer(res)
    A = A_ar(ϕ)
    plot(f_unitless, psd_vk ./ maximum(psd_vk), xscale=:log10, yscale=:log10)
    psd_fit = psd_ar.(f_unitless, Ref(ϕ))
    plot!(f_unitless, psd_fit ./ maximum(psd_fit))
end

A = A_ar(hcat(ol[2:end-1], ol[1:end-2]) \ ol[3:end])

ol = von_karman_turbulence(1000, 10)
begin
    state = []
    predictions_history = [0.0]
    for (i, f) in enumerate(ol)
        pushfirst!(state, f)
        if length(state) == 3
            deleteat!(state, 3)
            prediction = f # (A*state)[1]
        else
            prediction = f
        end
        if i < length(ol)
            push!(predictions_history, prediction)
        end
    end
    rms(predictions_history .- ol)
end

states = rand(2)
for i in 1:100
    push!(states, (A*states[end-1:end])[1])
end
plot(states)


# see if AR1 blows up??
begin
    x = rand()
    highest_so_far = x
    for _ in 1:10000
        x = 0.995 * x + rand(Normal(0, 1e-2))
        highest_so_far = max(highest_so_far, x)
    end
    highest_so_far
end