import pandas as pd
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

data = pd.read_csv('../result/task2/full_trajectory_mt-272.60_angle-66.00.csv', sep=' ')[::300]

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(-63051822463.82 - 5e8, -63051822463.82 + 5e8)
ax.set_ylim(141005424242.41 - 5e8, 141005424242.41 + 5e8)
line_satellite, = ax.plot([], [], 'b-', label='Траектория спутника')  # Синия линия для спутника
line_rocket, = ax.plot([], [], 'r-', label='Траектория ракеты')  # Красная линия для ракеты
point_satellite, = ax.plot([], [], 'bo')  # Точка для спутника
point_rocket, = ax.plot([], [], 'ro')  # Точка для ракеты


def init():
    line_satellite.set_data([], [])
    line_rocket.set_data([], [])
    point_satellite.set_data([], [])
    point_rocket.set_data([], [])
    return line_satellite, line_rocket, point_satellite, point_rocket


def update(frame):
    # Обновление данных для спутника
    x_data_satellite = data.x3[:frame]
    y_data_satellite = data.y3[:frame]
    point_satellite.set_data([data.x3.iloc[frame]], [data.y3.iloc[frame]])
    line_satellite.set_data(x_data_satellite, y_data_satellite)

    # Обновление данных для ракеты
    x_data_rocket = data.x[:frame]
    y_data_rocket = data.y[:frame]
    point_rocket.set_data([data.x.iloc[frame]], [data.y.iloc[frame]])
    line_rocket.set_data(x_data_rocket, y_data_rocket)

    return line_satellite, line_rocket, point_satellite, point_rocket


# Создание анимации
ani = FuncAnimation(fig, update, frames=len(data), init_func=init, blit=True, interval=0.00001)
ani.save('satellite_rocket_trajectory.mp4', writer='ffmpeg', fps=30)
# plt.legend()
# plt.show()
