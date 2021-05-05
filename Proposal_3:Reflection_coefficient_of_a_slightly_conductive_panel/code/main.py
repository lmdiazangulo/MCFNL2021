from mesh import Mesh
from solver import FDTD, Utilities, Source, Panel
from viewer import Animator
import copy
import numpy as np
from matplotlib import pyplot as plt

#Data input
sigma=0.04
permittivity=4
ddx=0.002
start_m=110
end_m=140
thickness=(end_m-start_m)*ddx
time=10e-8


malla1=Mesh(200,ddx,permittivity,sigma,start_m,end_m)
malla2=Mesh(200,ddx,1,0,start_m,end_m)
pulso=Source('gauss',40,12,20)

et1k1, et1k2, ex_film = FDTD(malla1,pulso,time).FDTDLoop(110,140)
e2tk1, e2tk2, _= FDTD(malla2,pulso,time).FDTDLoop(110,140)

Animator().animationex(ex_film,malla1)

#FFT results  
r, t = Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)
freq = Utilities().frequency(time,et1k1)

#Analytical Result
omega = np.linspace(1, 3e12, 100001) * 2 *np.pi
x = np.abs(Panel(thickness,  permittivity,   sigma).R(omega))
y = np.abs(Panel(thickness,  permittivity,   sigma).T(omega))

#Plot results
plt.figure()
plt.xlim((1,7e10))
plt.ylim((0,1))
plt.plot(omega, y, label='Exact T')
plt.plot(freq, t, label='FDTD T')
plt.plot(omega, x, label='Exact R')
plt.plot(freq, r, label='FDTD R')
plt.plot(freq, r*r+t*t, label=r'$R^2 +T^2$')
plt.grid()
plt.legend(loc='best', ncol=4)
plt.show()


#Animator().fftgraph(freq,r,t)
#Animator().fftexact(omega, x, y)