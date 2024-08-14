# %%
import numpy as np
from multiwfs.utils import rms
from multiwfs.dare import solve_dare
from multiwfs.dynamics import StateSpaceDynamics, StateSpaceObservation, simulate
from multiwfs.controller import Openloop, Integrator, LQG, MPC

np.random.seed(1)

dynamics = StateSpaceDynamics(
    np.array([[0.995, 0.0, 0.0], [0.0, 0.995, 0.0], [0.0, 0.0, 0.0]]), # state to state
    np.array([[0.0], [0.0], [1.0]]), # state to input
    np.array([[1e-2, 0.0, 0.0], [0.0, 1e-2, 0.0], [0.0, 0.0, 0.0]]), # state covariance
)

openloop = Openloop(p=dynamics.input_size)
integrator = Integrator(s=1, p=1, gain=0.3, leak=0.999)

observation = StateSpaceObservation(
    np.array([[1.0, 0.0, -1.0], [1.0, 1.0, -1.0]]), # state to measure
    np.array([[0.0], [0.0]]), # input to measure
    np.array([[1e-2, 0], [0, 1e-2]]) # measure covariance
)

# %%
# put this into the LQG class
lqg = LQG(dynamics, observation)
# Cost matrix to minimize error on the science WFS without caring about the auxiliary WFS
Ccost = np.copy(lqg.C)
Ccost[1] = 0.0
lqg.Q = Ccost.T @ Ccost
lqg.Pcon = solve_dare(dynamics.A, dynamics.B, lqg.Q, lqg.R, S=lqg.S)
lqg.L = -np.linalg.pinv(lqg.R + dynamics.B.T @ lqg.Pcon @ dynamics.B) @ (lqg.S.T + dynamics.B.T @ lqg.Pcon @ dynamics.A)
mpc = MPC(dynamics, observation, Q=lqg.Q, R=lqg.R)
sim = simulate(dynamics, observation, [openloop, integrator, lqg, mpc], nsteps=100);
# %%
rms(sim["LQG"]["noiseless_measurements"][:,0]) # error on the science WFS
# %%
rms(sim["LQG"]["noiseless_measurements"][:,1])
# error on the auxiliary WFS
# %%
rms(sim["MPC"]["noiseless_measurements"][:,0]) # error on the science WFS
# %%
rms(sim["MPC"]["noiseless_measurements"][:,1])

# %%
