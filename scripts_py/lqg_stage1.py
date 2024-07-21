# %%
import numpy as np
from multiwfs.dynamics import StateSpaceDynamics
from multiwfs.controller import Openloop, Integrator, LQG

np.random.seed(1)

dynamics = StateSpaceDynamics(
    np.array([[0.995]]), # state to state
    np.array([[0.0]]), # state to input
    np.array([[1.0], [1.0]]), # state to measure
    np.array([[1.0], [1.0]]), # input to measure
    np.array([[1e-2]]), # state covariance
    np.array([[1e-4, 0], [0, 1e-4]]) # measure covariance
)

openloop = Openloop(p=dynamics.input_size)
integrator = Integrator(s=dynamics.state_size, p=dynamics.input_size)
lqg = LQG(dynamics)
controllers = [openloop, integrator, lqg]
dynamics.simulate(controllers);
# %%
