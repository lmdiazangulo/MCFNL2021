import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import speed_of_light, epsilon_0, mu_0

eta_0 = np.sqrt(mu_0/epsilon_0)

class Panel: 
    def __init__(self, thickness, epsilon_r = 1.0, sigma = 0.0, mu_r = 1.0):
        self.thickness = thickness
        self.epsilon_r = epsilon_r
        self.mu_r = mu_r
        self.sigma = sigma

    def epsilon_c(self, omega):
        return self.epsilon_r*epsilon_0 - complex(0,1)*self.sigma/omega

    def mu_c(self, omega):
        return self.mu_r * mu_0

    def gamma(self, omega):
        return complex(0,1) * omega * \
            np.sqrt(self.epsilon_c(omega) * self.mu_c(omega))

    def eta(self, omega):
        return np.sqrt(self.mu_c(omega) / self.epsilon_c(omega))

    def phi(self, omega):
        gd  = self.gamma(omega) * self.thickness
        eta = self.eta(omega)
        return np.array([[np.cosh(gd),      np.sinh(gd) * eta], \
                         [np.sinh(gd) /eta, np.cosh(gd)      ]])

    def _den(self, omega):
        phi = self.phi(omega)
        return phi[0,0]*eta_0 + phi[0,1] + phi[1,0]*eta_0**2 + phi[1,1]*eta_0
        
    def T(self, omega):
        return  2*eta_0 / self._den(omega)

    def R(self, omega): 
        phi = self.phi(omega)
        return \
            (phi[0,0]*eta_0 + phi[0,1] - phi[1,0]*eta_0**2 - phi[1,1]*eta_0) / \
            self._den(omega)

plt.figure()
omega = np.linspace(1e2, 1e10, 1e2+1) * 2 * np.pi
plt.plot(omega, np.abs(Panel(10e-3,  5,   .0).R(omega)))
plt.plot(omega, np.abs(Panel(10e-3,  5,   .0).T(omega)))

plt.figure()
omega = np.logspace(2, 10, 1e2+1) * 2 * np.pi
plt.loglog(omega, np.abs(Panel(10e-6,  1, 50e3).T(omega)))
plt.loglog(omega, np.abs(Panel(100e-6, 1, 50e3).T(omega)))

plt.show()

print('=== Program finished ===')