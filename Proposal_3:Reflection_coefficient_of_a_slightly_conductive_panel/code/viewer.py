import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import rc

class Animator:

    def animationex(self, exanimation, malla):
        self.exanimation=exanimation
        self.malla=malla

        cb=malla.material()[1]
       
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set(xlim=(0, 200), ylim=(-1.2, 1.2))

        x = np.linspace(0, 200, 201)

        line = ax.plot(x, exanimation[0, :], color='k', lw=2)[0]                

        def animate(i):
            line.set_ydata(exanimation[i, :])

        plt.ylabel('E$_x$', fontsize='14')
       
        plt.plot((0.5 / cb - 1) / 3, 'k--',
                 linewidth=0.75) # The math on cb is just for scaling

        plt.text(170, 0.5, 'Eps = {}'.format(malla.epsilon_r),
                horizontalalignment='center')
        plt.text(170, -0.5, 'Cond = {}'.format(malla.sigma),
                horizontalalignment='center')
        plt.xlabel('FDTD cells')

        plt.subplots_adjust(bottom=0.25, hspace=0.45)
    

        anim=FuncAnimation(fig, animate, interval=1, frames=1500)
        
        plt.draw()
        plt.show()    

    def fftgraph(self, freq, r, t):
        self.freq=freq
        self.r=r
        self.t=t
        
        plt.plot(freq,r, label='R')
        plt.plot(freq,t, label='T')
        plt.plot(freq,r*r + (t*t), label='$R^2 + T^2 $')

        plt.xlim(0, 5e10)
        plt.ylim(-0.2, 1.2)
        plt.xlabel('Frequency ($ \omega $) ')
        
        plt.title('Reflectance and Transmittance')

        plt.legend()
        plt.show()
    
    def fftexact(self, omega, r, t):
        self.omega=omega
        self.r=r
        self.t=t
        
        plt.plot(omega,r, label='R')
        plt.plot(omega,t, label='T')
        plt.plot(omega,r*r + (t*t), label='$R^2 + T^2 $')

        plt.xlim(0, 5e10)
        plt.ylim(-0.2, 1.2)
        plt.xlabel('Frequency ($ \omega $) ')
        
        plt.title('Analytical Reflectance and Transmittance')

        plt.legend()
        plt.show()