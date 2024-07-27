"""
Various control strategies.
"""

from abc import ABC
from functools import partial
import numpy as np
import cvxpy as cp

from .dare import solve_dare

class Controller(ABC):
    def reset(self):
        pass

    def __call__(self, measurement):
        self.observe_law(measurement)
        return self.control_law()

    def observe_law(self, measurement):
        pass

    def control_law(self):
        pass

class Openloop(Controller):
    def __init__(self, p=2):
        self.name = "openloop"
        self.u = np.zeros((p,))
    
    def control_law(self):
        return self.u

class Integrator(Controller):
    def __init__(self, s, p, gain=0.5, leak=0.999):
        self.s = s
        self.p = p
        self.gain = gain
        self.leak = leak
        self.u = np.zeros((p,))
        self.state = np.zeros((s,))
        self.name = "integrator"

    def reset(self):
        self.u = np.zeros((self.p,))

    def observe_law(self, measurement):
        self.state = measurement[:self.p]
        
    def control_law(self):
        self.u = self.gain * self.state + self.leak * self.u
        return self.u

class LQG(Controller):
    def __init__(self, dyn, obs, name="LQG"):
        self.name = name
        self.A, self.B, self.C, self.D = dyn.A, dyn.B, obs.C, obs.D
        Q = obs.C.T @ obs.C
        R = obs.D.T @ obs.D
        S = obs.C.T @ obs.D
        self.x = np.zeros((dyn.state_size,))
        self.u = np.zeros((dyn.input_size,))
        self.Pobs = solve_dare(dyn.A.T, obs.C.T, dyn.W, obs.V)
        self.Pcon = solve_dare(dyn.A, dyn.B, Q, R, S=S)
        self.K = self.Pobs @ obs.C.T @ np.linalg.pinv(obs.C @ self.Pobs @ obs.C.T + obs.V)
        self.L = -np.linalg.pinv(R + dyn.B.T @ self.Pcon @ dyn.B) @ (S.T + dyn.B.T @ self.Pcon @ dyn.A)
        
    def measure(self):
        return self.C @ self.x + self.D @ self.u
        
    def predict(self):
        self.x = self.A @ self.x + self.B @ self.u

    def update(self, y):
        self.x = self.x + self.K @ (y - self.measure())
        
    def control_law(self):
        self.u = self.L @ self.x
        return self.u

    def observe_law(self, measurement):
        self.predict()
        self.update(measurement)
        
class MPC(Controller):
    def __init__(self, dyn, obs, name="MPC", horizon=1):
        self.name = name
        self.horizon = horizon
        self.A, self.B, self.C, self.D = dyn.A, dyn.B, obs.C, obs.D
        self.Q = obs.C.T @ obs.C
        self.R = obs.D.T @ obs.D
        self.S = obs.C.T @ obs.D
        self.state_size, self.input_size = dyn.state_size, dyn.input_size
        self.x = np.zeros((dyn.state_size,))
        self.x_opt = cp.Variable((dyn.state_size, horizon + 1))
        self.u = np.zeros((dyn.input_size,))
        self.u_opt = cp.Variable((dyn.input_size, horizon))
        self.Pobs = solve_dare(dyn.A.T, obs.C.T, dyn.W, obs.V)
        self.K = self.Pobs @ obs.C.T @ np.linalg.pinv(obs.C @ self.Pobs @ obs.C.T + obs.V)
        cost = 0
        constr = []
        self.x_curr = cp.Parameter((self.state_size,))
        self.x_curr.value = self.x
        for t in range(self.horizon):
            xt, xtp1, ut = self.x_opt[:,t], self.x_opt[:,t+1], self.u_opt[:,t]
            cost += cp.quad_form(xtp1, self.Q) + cp.quad_form(ut, self.R)
            constr += [xtp1 == self.A @ xt + self.B @ ut]
        constr += [self.x_opt[:, 0] == self.x_curr]
        self.problem = cp.Problem(cp.Minimize(cost), constr)
        
    def measure(self):
        return self.C @ self.x + self.D @ self.u
        
    def predict(self):
        self.x = self.A @ self.x + self.B @ self.u

    def update(self, y):
        self.x = self.x + self.K @ (y - self.measure())
        
    def control_law(self):
        # at the current time, you're constrained on state = the KF-recovered state
        # due to the separation principle, this is the best we can do
        # after that, we assume x[next] = Ax[curr] + Bu[curr]
        # LQG doesn't look at W or V anyway so that's fine
        # and we try and minimize x^T Qx + x^T Su + u^T Ru
        self.x_curr.value = self.x
        self.problem.solve()
        self.u = self.u_opt.value[:,0]
        return self.u

    def observe_law(self, measurement):
        self.predict()
        self.update(measurement)