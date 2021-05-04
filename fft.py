import numpy as np
import matplotlib.pyplot as plt

s0=10e-6
t0=10*s0
num=1000+1

t = np.linspace(0,500e-6,num)
gauss=np.exp(-np.power(t-t0,2)/2/s0**2)
# np.power se usa porque tengo un array no un numero
# plt.figure('Dominio del tiempo')
# plt.plot(t,gauss)
# plt.show()

plt.figure('Dominio de la frecuencia')
n=gauss.size
f = np.fft.fftfreq(n)
gaussf=np.fft.fft(gauss)
plt.plot(f,np.absolute(gaussf))
plt.show()