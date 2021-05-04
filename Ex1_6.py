import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0, mu_0

# Ejercicio 1.6 Coeficientes de transmision y reflexion

eta_0=np.sqrt(mu_0/epsilon_0)

class Panel:
    def __init__(self, thickness, epsilon_r, sigma, mu_r):
        self.thickness = thickness
        self.epsilon_r = epsilon_r
        self.mu_r = mu_r
        self.sigma = sigma
    
    def epsilon_c(self, omega):
        return self.epsilon_r*epsilon_0 - complex(0,1)*self.sigma/omega

    def mu_c(self, omega):
        return self.mu_r * mu_0

    def gamma(self, omega):
        return complex(0,1) * omega * np.sqrt(self.mu_c(omega) * self.epsilon_c(omega))
    
    def eta(self, omega):
        return np.sqrt(self.mu_c(omega)/self.epsilon_c(omega))

    def phi(self, omega):
        gd=self.gamma(omega) * self.thickness
        eta=self.eta(omega)

        return np.array([ [np.cosh(gd), np.sinh(gd) * eta], \
            [np.sinh(gd) / eta, np.cosh(gd)]  ])

    def _den(self,omega):
        phi=self.phi(omega)
        return phi[0,0]*eta_0 + phi[0,1] + phi[1,0]*eta_0*eta_0 + phi[1,1]*eta_0

    def T(self,omega):
        return 2*eta_0 / self._den(omega)

    def R(self,omega):
        phi=self.phi(omega)
        return (phi[0,0]*eta_0 + phi[0,1] - phi[1,0]*eta_0*eta_0 - phi[1,1]*eta_0) / self._den(omega)
    


Apartado1=Panel(100e-6,5,0,1)
Apartado2=Panel(10e-3,5,0,1)
Apartado3=Panel(100e-6,1,10e3,1)
Apartado4=Panel(10e-3,1,10e3,1)


omega=np.linspace(1e3,1e9,1000)

plt.figure('Apartado 1')
plt.subplot(2,1,1)
plt.title('Coeficiente de reflexion')
plt.plot(2 * np.pi * omega, np.abs(Apartado1.R(omega)))
plt.subplot(2,1,2)
plt.title('Coeficiente de transmision')
plt.plot(2 * np.pi * omega, np.abs(Apartado1.T(omega)))


plt.figure('Apartado 2')
plt.subplot(2,1,1)
plt.title('Coeficiente de reflexion')
plt.plot(2 * np.pi * omega, np.abs(Apartado2.R(omega)))
plt.subplot(2,1,2)
plt.title('Coeficiente de transmision')
plt.plot(2 * np.pi * omega, np.abs(Apartado2.T(omega)))


plt.figure('Apartado 3')
plt.subplot(2,1,1)
plt.title('Coeficiente de reflexion')
plt.plot(2 * np.pi * omega, np.abs(Apartado3.R(omega)))
plt.subplot(2,1,2)
plt.title('Coeficiente de transmision')
plt.plot(2 * np.pi * omega, np.abs(Apartado3.T(omega)))


plt.figure('Apartado 4')
plt.subplot(2,1,1)
plt.title('Coeficiente de reflexion')
plt.plot(2 * np.pi * omega, np.abs(Apartado4.R(omega)))
plt.subplot(2,1,2)
plt.title('Coeficiente de transmision')
plt.plot(2 * np.pi * omega, np.abs(Apartado4.T(omega)))
plt.show()


    


