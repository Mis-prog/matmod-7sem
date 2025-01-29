from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import sys

b = 0

data = np.genfromtxt(f'../result/symplectic.txt', delimiter=' ',max_rows=700)
print(len(data))
fig = plt.figure()
my = max([max(dat) for dat in data])
my = my / 2
ax = plt.axes(xlim=(0, len(data[0])), ylim=(-4, 4))

x = list(range(len(data[0])))
line, = ax.plot(x, data[0], lw=2)
i = 1
st = 0
fr = 1000


def init():
    line.set_data(x, data[0])
    return line,


def animate(i):
    try:
        line.set_data(x, data[st + i])
        i += 1
        return line,
    except:
        sys.exit(1)


anim = FuncAnimation(fig, animate,
                     frames=fr, interval=20, blit=True)
# plt.show()

anim.save('3_fpu_betta.mp4', writer='ffmpeg', fps=25)