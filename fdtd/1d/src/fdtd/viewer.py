import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Animator:

    def __init__(self, mesh, probe):
        
        ids = probe["indices"]
        gridE = mesh.pos[ids[0]:ids[1]]

        probeTime = probe["time"][:]
        values    = probe["values"][:]


        fig = plt.figure(figsize=(8,4))
        ax1 = fig.add_subplot(1, 1, 1)
        ax1 = plt.axes(xlim=(gridE[0], gridE[-1]), ylim=(-1.1, 1.1))
        ax1.grid(color='gray', linestyle='--', linewidth=.2)
        ax1.set_xlabel('X coordinate [m]')
        ax1.set_ylabel('Field')
        line1,    = ax1.plot([], [], '.', markersize=1)
        timeText1 = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)

        def init():
            line1.set_data([], [])
            timeText1.set_text('')
            return line1, timeText1

        def animate(i):
                line1.set_data(gridE, values[i][:])
                timeText1.set_text('Time = %2.1f [ns]' % (probeTime[i]*1e9))
                return line1, timeText1

        _ = animation.FuncAnimation(fig, animate, init_func=init,
            frames=len(probeTime), interval=10, blit=True)

        plt.show()