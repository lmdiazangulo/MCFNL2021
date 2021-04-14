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
    return A*np.exp(-((x-c*t)/(sig*np.sqrt(2)))**2)
    
tEnd=100
xEnd=100
c=1
A=1
sig=1

t=np.linspace(0,tEnd,50)
x=np.linspace(0,xEnd,5000)

fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(1, 2, 1)
ax1 = plt.axes(xlim=(x[0], x[-1]), ylim=(-1.1, 1.1))
line1,    = ax1.plot([], [], 'o', markersize=1)
timeText1 = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)

def init():
    line1.set_data([], [])
    timeText1.set_text('')
    return line1, timeText1

def animate(i):
    line1.set_data(x, Gaussian(t[i], x, A, sig, c))
    timeText1.set_text('Time = %2.1f' % (t[i]))
    return line1, timeText1

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(t), interval=50, blit=True)

plt.show()


