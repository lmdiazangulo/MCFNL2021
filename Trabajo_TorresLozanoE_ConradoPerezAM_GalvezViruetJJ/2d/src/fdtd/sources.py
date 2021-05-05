import numpy as np
import scipy.constants as sp

def TMn_Ex(t, x, y, n, intens, freq, mu, epsilon, d):
    '''
    Gives the Ex field of the TMn mode of a parallel plate waveguide.
    '''
    kc = n*np.pi/d
    beta = np.sqrt((freq**2*mu*epsilon - kc**2).astype(np.complex))
    return np.real(intens\
        * np.exp(- 1j* beta * x) * np.sin(kc * y) \
        * np.exp(1j *freq * t)) 

def TMn_Ey(t, x, y, n, intens, freq, mu, epsilon, d):
    '''
    Gives the Ey field of the TMn mode of a parallel plate waveguide.
    '''
    kc = n*np.pi/d
    beta = np.sqrt((freq**2*mu*epsilon - kc**2).astype(np.complex))
    return np.real(- intens * 1j * beta / kc \
        * np.exp(- 1j* beta * x) * np.cos(kc * y) \
        * np.exp(1j *freq * t))

def TMn_Hz(t, x, y, n, intens, freq, mu, epsilon, d):
    '''
    Gives the Hz field of the TMn mode of a parallel plate waveguide.
    '''
    kc = n*np.pi/d
    beta = np.sqrt((freq**2*mu*epsilon - kc**2).astype(np.complex))
    return np.real(intens * 1j * freq * epsilon / kc \
        * np.exp(- 1j* beta * x) * np.cos(kc * y) \
        * np.exp(1j *freq * t))

def gaussian(x, delay, spread):
    '''
    Gives the value of a gaussian pulse.
    '''
    return np.exp(-((x-delay)**2 / (2*spread**2)))

def step(x, xlim):
    '''
    Gives the value x of the Heaviside function with the step located in xlim.
    '''
    return x<xlim

def movingGaussian(x,y,t,c,center,A,spread):
    '''
    Gives the value of a moving gaussian function.
    '''
    return A*np.exp(-(((x-center)-c*t)**2 /(2*spread**2)))

def fc_f(n, d, mu, epsilon, freq, _has_run=[]):
    '''
    Calculates and shows the value of the cut off frequency of the TM mode 
    just once
    '''
    if _has_run: return
    fc = n / (2 * d * np.sqrt(mu * sp.mu_0 * epsilon * sp.epsilon_0))
    print(f'Higher cut off frequency of the TM mode in the selected materials = {np.amax(fc):e}')

    if (freq < fc).any(): print("Warning: Selected frequency is below cut off frequency in some material selected")
    _has_run.append(1)