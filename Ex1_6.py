import numpy as np
import matplotlib.pyplot as plt
import math as m
import cmath as cm
# Ejercicio 1.6 Coeficientes de transmision y reflexion

def eps(epsr,sigma,f):
    ''' Funcion que devuelve la permitividad efectiva complega
    con parametros la perm. relativa, conductividad y frec. '''
    eps0=8.85e-12 # C^2/(N*m^2)
    perm=epsr*eps0-sigma/(2*m.pi*f)*1j
    return perm

# Constantes

eps0=8.85e-12 # C^2/(N*m^2)
mu0=4*m.pi*1e-7 # T*m*A^-1

# Propiedades del medio 1
epsr=1
d=10e-3
sigma=10e3

f=np.linspace(1e3,1e9,1000000)
T=np.zeros(len(f),dtype=np.complex_)
R=T

for i in [0,1,len(f)-1]:

    perm=eps(epsr,sigma,f[i])
    gamma=1j*(2*m.pi)*f[i]*cm.sqrt(mu0*perm)
    eta=cm.sqrt(mu0/perm)

    phi00=cm.cosh(gamma*d)
    phi11=phi00
    phi01=eta*cm.sinh(gamma*d)
    phi10=1/eta*cm.sinh(gamma*d)
    
    T[i]=2*eta/(phi00*eta+phi01+phi10*eta*eta+phi11*eta)
    R[i]=(phi00*eta+phi01-phi10*eta*eta-phi11*eta)/(phi00*eta+phi01+phi10*eta*eta+phi11*eta)

fig, (ax1, ax2) = plt.subplots(2)
ax1.set_title('Coeficiente transmision')
ax1.plot(f,T)
ax2.set_title('Coeficiente reflexion')
ax2.plot(f,R)
plt.show()



    


