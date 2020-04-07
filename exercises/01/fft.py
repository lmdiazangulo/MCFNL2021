import numpy as np
import math
import matplotlib.pyplot as plt

def gaussian(t, delay, spread):
    return np.exp(- np.power(t-delay, 2) / (2.0 * np.power(spread, 2)) )


s0 = 10e-6
t0 = 10 * s0
N  = int (1e2+1)

tIni = 0
tEnd = 500e-6

t = np.linspace(tIni, tEnd, num = N, endpoint=True)
f = gaussian(t, t0, s0)

fq = np.fft.fftfreq(N) / (t[1]-t[0])
f_fq = np.fft.fft(f)

plt.figure()
plt.plot(t, f)

plt.figure()
plt.plot(np.fft.fftshift(fq), np.fft.fftshift(np.abs(f_fq)))

plt.show()

print('=== Program finished ===')