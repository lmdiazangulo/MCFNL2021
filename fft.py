import numpy as np
import matplotlib.pyplot as plt

s0 = 10e-6
t0 = 10*s0

t = np.linspace(0, 500e-6, num=(1000 + 1))
gauss = np.exp( - np.power(t - t0, 2) / 2 / s0**2)

plt.plot(t, gauss)
plt.show()
