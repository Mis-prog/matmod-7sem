import pandas as pd
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

data = pd.read_csv('../result/task2/trajectory_mt-2.05263_angle-1.csv', sep=' ')[::200]

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(-50368219856.43 - 2e7, -50368219856.43 + 2e7)
ax.set_ylim(-219503615669.65 - 2e7, -219503615669.65+ 2e7)
line_satellite, = ax.plot([], [], 'b-', label='Траектория спутника')  # Синия линия для спутника
line_rocket, = ax.plot([], [], 'r-', label='Траектория ракеты')  # Красная линия для ракеты
point_satellite, = ax.plot([], [], 'bo')  # Точка для спутника
point_rocket, = ax.plot([], [], 'ro')  # Точка для ракеты
point_planeta = ax.scatter(-50368219856.43,-219503615669.65,color='orange',label='Марс')

def init():
    
    line_satellite.set_data([], [])
    point_satellite.set_data([], [])
    return line_satellite, point_satellite


def update(frame):
    # Обновление данных для спутника
    x_data_satellite = data.x3[:frame]
    y_data_satellite = data.y3[:frame]
    point_satellite.set_data([data.x3.iloc[frame]], [data.y3.iloc[frame]])
    line_satellite.set_data(x_data_satellite, y_data_satellite)

    # Обновление данных для ракеты
    # x_data_rocket = data.x[:frame]
    # y_data_rocket = data.y[:frame]
    # point_rocket.set_data([data.x.iloc[frame]], [data.y.iloc[frame]])
    # line_rocket.set_data(x_data_rocket, y_data_rocket)

    return line_satellite,  point_satellite


# Создание анимации
ani = FuncAnimation(fig, update, frames=len(data), init_func=init, blit=True, interval=1)
ani.save('satellite_rocket_trajectory.mp4', writer='ffmpeg', fps=30)
# plt.grid()
# plt.legend()
# plt.show()
