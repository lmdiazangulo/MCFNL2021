#fft.py
import numpy as np
import scipy.stats as ss
import math
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy import optimize




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
freq = np.fft.fftfreq(t.shape[-1])
gft=np.fft.fft(gaussian(t,t0,s0))
plt.subplot(1, 2, 2)
plt.plot(freq,abs(gft),color='r')
plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.ylabel('G(w) ')
plt.xlabel('w')
plt.show()
