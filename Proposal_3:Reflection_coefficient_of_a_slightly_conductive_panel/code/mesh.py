import numpy as np
import scipy.constants as sp
from math import pi, sin, exp

class Mesh:
    def __init__(self, ncells, ddx, epsilon_r, sigma, start_m, end_m):    
        self.ncells=ncells
        self.ddx=ddx
        self.epsilon_r= epsilon_r
        self.sigma= sigma
        self.start_m= start_m
        self.end_m= end_m
    
    def material(self):
        
        dt=self.ddx / (2*sp.c)

        ca = np.ones(self.ncells+1)
        cb = np.ones(self.ncells+1) * 0.5
        
        eaf = dt * self.sigma / (2 * sp.epsilon_0 * self.epsilon_r)
        ca[self.start_m : self.end_m] = (1 - eaf ) / (1 + eaf )
        cb[self.start_m : self.end_m] = 0.5 / (self.epsilon_r * (1 + eaf ))

        return ca, cb
