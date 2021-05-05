import numpy as np
import copy
from math import pi, sin, exp
import scipy.constants as sp


class FDTD:
    def __init__(self, mesh, pulse, time):
        self.mesh=mesh
        self.pulse=pulse
        self.time=time

    def boundarymur(self, ex, ex_old): 
        ncells, ddx= self.mesh.ncells, self.mesh.ddx
        dt=self.mesh.ddx / (2*sp.c)

        c_bound=(sp.c*dt-ddx)/(sp.c*dt+ddx)

        ex[0]=ex_old[1] + c_bound * (ex[1]-ex_old[0])
        ex[ncells]=ex_old[ncells-1] + c_bound * (ex[ncells-1]-ex_old[ncells])


    def FDTDLoop(self,k1,k2):
        
        dt=self.mesh.ddx / (2*sp.c)
        nsteps= int(self.time  / dt)

        # COMENTAR: Mejor quitar nsteps, no guardar siempre todo...
        ex=np.zeros(self.mesh.ncells+1)
        hy=np.zeros(self.mesh.ncells+1)
        ex_old=np.zeros(self.mesh.ncells+1)

        ex_save_k1=np.empty(nsteps+1)
        ex_save_k2=np.empty(nsteps+1)

        ex_save_film=np.empty((nsteps+1,self.mesh.ncells+1))
        
        ca=self.mesh.material()[0][1:-1]
        cb=self.mesh.material()[1][1:-1]
       
        for time_step in range(1, nsteps + 1):
            ex_old=copy.deepcopy(ex)

            ex[1:-1] = ca * ex[1:-1] + cb * (hy[:-2] - hy[1:-1])
            
            #Guardo los valores a representar
            ex_save_film[time_step][:]=ex[:]
            
            #Guardo los valores para calcular la transformada
            ex_save_k1[time_step]=ex[k1]
            ex_save_k2[time_step]=ex[k2]
           
            ex[self.pulse.k_ini] +=  0.5*self.pulse.pulse(time_step) 
            
            self.boundarymur(ex,ex_old) 
            
            
            hy[:-1] = hy[:-1] + 0.5 * (ex[:-1] - ex[1:])   

            t= time_step+1/2
            hy[self.pulse.k_ini] += 0.25* self.pulse.pulse(t) 
            hy[self.pulse.k_ini-1] += 0.25* self.pulse.pulse(t)   

                       

        return ex_save_k1, ex_save_k2,  ex_save_film




class Source:
    def __init__(self, sourcetype, t_0, s_0, k_ini):
        self.sourcetype=sourcetype
        self.t_0=t_0
        self.s_0=s_0
        self.k_ini=k_ini

    def pulse(self, time):
        
        self.time=time
        
        if self.sourcetype == 'gauss':
            pulse = exp(-0.5*( (self.t_0 - time) / self.s_0 )**2)
        
        return pulse




#Clase para la Trasformada Rápida de Fourier
# COMENTAR: Esto es mas un namespace que una clase. 
# COMENTAR: Cuanto menos estado, mejor
class Utilities:

    def FFT(self,e1tk1_total,e2tk1,e1tk2,e2tk2):
        
        #Hay que cancelar la parte incidente
        e1tk1_reflected = e1tk1_total - e2tk1  
        
        e1wk1=np.fft.fft(e1tk1_reflected)
        e2wk1=np.fft.fft(e2tk1)

        e1wk2=np.fft.fft(e1tk2)
        e2wk2=np.fft.fft(e2tk2)
    
        R=np.abs(e1wk1) / np.abs(e2wk1)
        T=np.abs(e1wk2) / np.abs(e2wk2)
        
        
        return  R, T
    

    def frequency(self,time,e1tk1):

        N=len(e1tk1)

        w= ((2.0*np.pi)/time) * np.arange(N)         

        return w

eta_0 = np.sqrt(sp.mu_0/sp.epsilon_0)
class Panel: 
    def __init__(self, thickness, epsilon_r = 1.0, sigma = 0.0, mu_r = 1.0):
        self.thickness = thickness
        self.epsilon_r = epsilon_r
        self.mu_r = mu_r
        self.sigma = sigma

    def epsilon_c(self, omega):
        return self.epsilon_r*sp.epsilon_0 - complex(0,1)*self.sigma/omega

    def mu_c(self, omega):
        return self.mu_r * sp.mu_0

    def gamma(self, omega):
        return complex(0,1) * omega * \
            np.sqrt(self.epsilon_c(omega) * self.mu_c(omega))

    def eta(self, omega):
        return np.sqrt(self.mu_c(omega) / self.epsilon_c(omega))

    def phi(self, omega):
        gd  = self.gamma(omega) * self.thickness
        eta = self.eta(omega)
        return np.array([[np.cosh(gd),      np.sinh(gd) * eta], \
                         [np.sinh(gd) /eta, np.cosh(gd)      ]])

    def _den(self, omega):
        phi = self.phi(omega)
        return phi[0,0]*eta_0 + phi[0,1] + phi[1,0]*eta_0**2 + phi[1,1]*eta_0
        
    def T(self, omega):
        return  2*eta_0 / self._den(omega)

    def R(self, omega): 
        phi = self.phi(omega)
        return \
            (phi[0,0]*eta_0 + phi[0,1] - phi[1,0]*eta_0**2 - phi[1,1]*eta_0) / \
            self._den(omega)

