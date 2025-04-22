using multiwfs
using NPZ
using CairoMakie

sim = simgen_ichpf(0.4, 1.4, 15.0, VonKarman(0.001, 0.25), 50.0; leak=0.995)

atm = npzread("data/olt/olt_atm.npy")
fastncp = npzread("data/olt/olt_fastncp.npy")
slowncp = npzread("data/olt/olt_slowncp.npy")
fastnoise = npzread("data/olt/olt_fastnoise.npy")
slownoise = npzread("data/olt/olt_slownoise.npy")

slowphase_ts_atm = timedomain_closedloop(sim, (atm=atm, fastncp=fastncp*0, slowncp=slowncp*0, fastnoise=fastnoise*0, slownoise=slownoise*0))
slowphase_ts_fastncp = timedomain_closedloop(sim, (atm=atm*0, fastncp=fastncp, slowncp=slowncp*0, fastnoise=fastnoise*0, slownoise=slownoise*0))
slowphase_ts_slowncp = timedomain_closedloop(sim, (atm=atm*0, fastncp=fastncp*0, slowncp=slowncp, fastnoise=fastnoise*0, slownoise=slownoise*0))

_, psd_atm = genpsd(atm, 1000.0)
_, psd_fastncp = genpsd(fastncp, 1000.0)
_, psd_slowncp = genpsd(slowncp, 1000.0)
freq, psd_cl_atm = genpsd(slowphase_ts, 1000.0)
freq, psd_cl_fastncp = genpsd(slowphase_ts_fastncp, 1000.0)
freq, psd_cl_slowncp = genpsd(slowphase_ts_slowncp, 1000.0)

begin
    fig = Figure()
    ax = Axis(fig[1,1], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)")
    lines!(ax, freq, psd_cl_atm ./ psd_atm)
    #lines!(ax, freq, psd_cl_fastncp ./ psd_fastncp)
    #lines!(ax, freq, psd_cl_slowncp ./ psd_slowncp)
    lines!(ax, sim.fr, abs2.(phi_to_X.(sim.sT, Ref(sim))), linestyle=:dash, color=:black)
    #lines!(ax, sim.fr, abs2.(Lfast_to_X.(sim.sT, Ref(sim))), linestyle=:dash, color=:black)
    #lines!(ax, sim.fr, abs2.(Lslow_to_X.(sim.sT, Ref(sim))), linestyle=:dash, color=:black)
    CairoMakie.ylims!(ax, 1e-5, 2e3)
    fig
end