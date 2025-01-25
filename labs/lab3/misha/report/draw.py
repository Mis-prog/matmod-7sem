from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

b = 0

data = np.genfromtxt(f'../result/speed_{b}.000000', delimiter=' ')
print(data)
fig = plt.figure()
my=max([max(dat) for dat in data])
my=my/2
ax = plt.axes(xlim=(0, len(data[0])), ylim=(-4, 4))

x=list(range(len(data[0])))
line, = ax.plot(x, data[0], lw=2)
i=1
st=0
fr=1000
def init():
    line.set_data(x, data[0])
    return line,
def animate(i):
    line.set_data(x, data[st+i])
    i+=1
    return line,


anim = FuncAnimation(fig, animate,
                    frames=fr, interval=50, blit=True)

title = f'soliton a = {0.7} b = {b}.gif'

plt.show()
