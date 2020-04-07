import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

x = np.linspace(-5, 5, num=101)
t = np.linspace(0, 10, num=101)

f = np.zeros((len(t), len(x))) 
f += 0.5

# --- Creates animation ---
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
    line1.set_data(x, f[i,:])
    timeText1.set_text('Time = %2.1f' % (t[i]*1e9))
    return line1, timeText1

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(t), interval=50, blit=True)

plt.show()
