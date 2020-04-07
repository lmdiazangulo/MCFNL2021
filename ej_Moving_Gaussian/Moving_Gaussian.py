# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 13:37:26 2020

@author: María Pedrosa Bustos (mpedrosab@gmail.com)
"""

# =============================================================================
# Gaussian function moving with speed c
# 
# =============================================================================

import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt


def Gaussian(t,x,A,sig,c):
    '''
        Computes Gaussian function
    
    '''
    f=A*np.exp(-((x-c*t)/(sig*np.sqrt(2)))**2)
    plt.plot()
    return f

def PlotMovie(t,x,A,sig,c) :
    f=Gaussian(t,x,A,sig,c)
    ax.clear()
    graf=ax.plot(x,f)
    return graf
    
tEnd=10
xEnd=100
c=1
A=1
sig=1

t=np.linspace(0,tEnd,50)
x=np.linspace(0,xEnd,5000)

#%% Plot
# Graphic 
fig,ax=plt.subplots(figsize=(12,7))
ax.set_xlabel("x")

# Movie properties
Writer = animation.writers['ffmpeg']
writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1,codec="libx264")

# Play movie
peli=animation.FuncAnimation(fig,PlotMovie,fargs=(x,A,sig,c),frames=t,interval=30,repeat=True)



