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
x_values =np.linspace(0,500e-6,num=(1000+1))
plt.plot(x_values,gaussian(x_values,t0,s0))


plt.show()