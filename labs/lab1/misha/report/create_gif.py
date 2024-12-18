import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data=pd.read_csv('../res_task2/full_trajectory.csv',sep=' ')
end = 100000
frame_step=1000
fig, ax = plt.subplots()

ax.set_xlim(-50368219856.43 - 1e12, -50368219856.43 + 1e12)
ax.set_ylim(-219503615669.65 - 1e12, -219503615669.65 + 1e12)

point_mars = ax.scatter(-50368219856.43, -219503615669.65, color='red', label='Марс')
line_rocket, = ax.plot([], [], color='orange', label='Ракета')
line_phobos, = ax.plot([], [], color='green', label='Фобос')
line_mars, = ax.plot([], [], color='blue', label='Марс')
ax.legend(loc='best')
ax.grid()

def update(frame):
    idx = frame * frame_step  # Расчёт текущего индекса
    line_rocket.set_data(data.x[:idx], data.y[:idx])
    line_phobos.set_data(data.x3[:idx], data.y3[:idx])
    line_mars.set_data(data.x2[:idx], data.y2[:idx])
    return line_rocket, line_phobos, line_mars


ani = FuncAnimation(
    fig, update, frames=range(0, end // frame_step), interval=50, blit=True
)

ani.save('animation.gif', writer='pillow')