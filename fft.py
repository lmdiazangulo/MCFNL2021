#fft.py
import numpy as np
import scipy.stats as ss
import math
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.fft import fft, fftfreq



s0=10e-6
t0=10*s0
def gaussian(x,mu,sig):
    return np.exp(-np.power(x-mu,2.)/(2*np.power(sig,2.)))
t =np.linspace(0,500e-6,num=(1000+1))
plt.subplot(1, 2, 1)
plt.plot(t,gaussian(t,t0,s0))
plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.ylabel('g(t) ')
plt.xlabel('t')
plt.title('Funcion Gaussiana')
plt.xlim(0,200e-6)
#---- transformada de Fourier-------

# Number of sample points
N = 1001

# sample spacing

T = 500e-6 / 1001.0

freq = fftfreq(N, T)[:N]

gft= fft(gaussian(t,t0,s0))

plt.subplot(1, 2, 2)
plt.plot(freq,abs(gft),color='r')
plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.ylabel('G(w) ')
plt.xlabel('w')
plt.show()
