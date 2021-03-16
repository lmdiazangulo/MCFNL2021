""" Un pequeño script para relizar la transformada rápida de Fourier """
import timeit
import numpy as np 
from scipy.fftpack import fft, fftfreq
import matplotlib.pyplot as plt

gaussian = lambda x,x0=0,σ=1: np.exp(-(x-x0)**2/(σ**2*2))/(np.sqrt(2*np.pi)*σ)
s0 = 10*1e-6 # μs 
fgauss = lambda x: gaussian(x,10*s0,s0)
fsin = lambda x, w = 1: np.sin(w*x)

num = (1000 + 1)
def Fourier_transform_plot(f,T,N = 100, k = 1, a = 0):
    """ Funcion que realiza la transformada de Fourier directa de la función f y plotea
    la función en el dominio del tiempo y la frecuencia.
    Inputs: 
        - f: Función en la que se quiere usar la transformada de Fourier directa.
        - T: Espaciado temporal, 1/2T frecuencia máxima.
        - N: N de puntos, el tiempo máximo es TN.
        - k: Parametro que permite expandir el intervalo de tiempo donde se 
        plotea la función dependiente del tiempo.
        - a: Tiempo inicial 
    Outputs:
        - Tupla cuyo primer elemento es el array resultado con números
    complejos y cuyo segundo es el módulo de estos números. """



    # Definimos el mallado temporal para transformada.
    t = np.linspace(a, N*T-a, N)

    # Mallado temporal para representación y imagines del mallado.
    t0 = list(map(lambda x, k0 = k: x*k0, t))
    f0 = [f(i) for i in t0]

    # Calculamos la transformada rápida de Fourier.
    fft0 = fft(f(t))
    xf = fftfreq(N, T)[:N//2]  # frecuencias positivas.

    # Calculamos el módulo del resultado.
    modfft = list(map(lambda x: 2/N * abs(x),fft0))

    # Ploteo de la función dependiente del tiempo.
    plt.plot(t0,f0)
    plt.show()

    # Plotep de la función dependiente de la frecuencia.
    plt.plot(xf, modfft[:N//2])
    plt.show()

    return (fft,modfft)



(_,fftgauss) = Fourier_transform_plot(fgauss,500e-6/1001,1001) 


(_,fftsin) = Fourier_transform_plot(lambda t: np.sin(0.1*2*np.pi*t),0.5,int(250),k = 1)
