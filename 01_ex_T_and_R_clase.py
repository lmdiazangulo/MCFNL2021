import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import speed_of_light, epsilon_0, mu_0

eta_0 = np.sqrt(mu_0/epsilon_0)

def _denominator(phi):
    return phi[0,0]*eta_0 + phi[0,1] + phi[1,0]*eta_0**2 + phi[1,1]*eta_0

fIni = 1e3
fEnd = 1000e9
N = 1e2 + 1

omega = np.linspace(fIni, fEnd, N) * 2 * np.pi

epsilon_r = 5.0
sigma     = 0.0
d         = 100e-6

epsilon_c = epsilon_r * epsilon_0 + sigma / (np.complex(0,1)*omega)
mu = 1.0 * mu_0
gamma = np.complex(0, 1) * omega * np.sqrt(mu * epsilon_c)

eta = np.sqrt(mu / epsilon_c) 
gd = gamma*d

phi = np.array([[np.cosh(gd), eta*np.sinh(gd)], \
                [np.sinh(gd)/eta, np.cosh(gd)]])

T = 2 * eta_0 / _denominator(phi)
R = (phi[0,0]*eta_0 + phi[0,1] - phi[1,0]*eta_0**2 - phi[1,1]*eta_0 ) / \
    _denominator(phi)


plt.plot(omega/2/np.pi, np.abs(T))
plt.plot(omega/2/np.pi, np.abs(R))
plt.show()

print('=== Program finished ===')