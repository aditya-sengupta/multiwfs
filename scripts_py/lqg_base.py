# %%
import numpy as np
from multiwfs import LQG

lqg = LQG(
    np.array([[0.995]]),
    np.array([[1.0]]),
    np.array([[1.0]]),
    np.array([[0.0]]),
    np.array([[1e-2]]),
    np.array([[1e-4]])
)

lqg.simulate();
# %%
