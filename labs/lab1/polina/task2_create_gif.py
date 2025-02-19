import pandas as pd
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

class Constants:
    G = 6.67e-11  # гравитационная постоянная
    M1 = 2.0e30   # масса звезды (кг)
    M2 = 5.7e26   # масса планеты (кг)
    M3 = 1.4e23   # масса спутника (кг)
    R1 = 696340e3 # радиус звезды (м)
    R2 = 60270e3   # радиус планеты (м)
    R3 = 2575e3   # радиус спутника (м)
    R12 = 1429e9   # начальное расстояние звезда-планета (м)
    R23 = 1222e6   # начальное расстояние планета-спутника (м)
    U2 = 9.7e3     # начальная скорость планеты (м/с)
    U3 = 5.57e3   # начальная скорость спутника (м/с)
    T = 4100.0    # время работы двигателя (с)
    H = 1000e3     # высота орбиты (м)
    M0 = 95.0     # масса полезной нагрузки (кг)
    U = 2500.0    # скорость истечения (м/с)
    koef = 0.001
    
r12y0, r12x0 = 0 , Constants.R1 + Constants.R12 + Constants.R2

data = pd.read_csv('trajectory_mt-5392.25708722_angle-260.7663678.csv', sep=' ')[::200]

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_xlim(r12x0- 2e9, r12x0 + 2e9)
ax.set_ylim(r12y0 - 2e9,r12y0+ 2e9)
line_satellite, = ax.plot([], [], 'b-', label='Траектория спутника')  # Синия линия для спутника
line_rocket, = ax.plot([], [], 'r-', label='Траектория ракеты')  # Красная линия для ракеты
point_satellite, = ax.plot([], [], 'bo')  # Точка для спутника
point_rocket, = ax.plot([], [], 'ro')  # Точка для ракеты
point_planeta = ax.scatter(r12x0, r12y0,color='orange',label='Сатурн')

def init():
    
    line_satellite.set_data([], [])
    point_satellite.set_data([], [])
    line_rocket.set_data([], [])
    point_rocket.set_data([], [])
    return line_satellite, point_satellite,line_rocket, point_rocket


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

    return line_satellite,  point_satellite,line_rocket , point_rocket


# Создание анимации
ani = FuncAnimation(fig, update, frames=len(data), init_func=init, blit=True, interval=1)
# plt.grid()
# plt.legend()
# plt.show()

ani.save('satellite_rocket_trajectory.mp4', writer='ffmpeg', fps=30)