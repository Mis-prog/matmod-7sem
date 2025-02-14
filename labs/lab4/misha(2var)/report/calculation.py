import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D

# Параметры системы
a, b = -0.5, 0.5

# Уравнения системы
def system(t, state):
    x, y, z = state
    dx = a * (y - x)
    dy = b * y - x * z
    dz = -3 * z + x * y
    return [dx, dy, dz]

static_points = np.array([
    [np.sqrt(3 * b), np.sqrt(3 * b), b],
    [-np.sqrt(3 * b), -np.sqrt(3 * b), b],
    [0, 0, 0]
])

# Начальные условия
# state0 = [1.0, 1.0, 1.0]
state0 = [0.1, 0.1, 0.1]
t_span = [0, 100]
t_eval = np.linspace(t_span[0], t_span[1], 10000)

# Решение системы
sol = solve_ivp(system, t_span, state0, t_eval=t_eval, method='RK45')

# Построение фазового портрета
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(static_points[:, 0], static_points[:, 1], static_points[:, 2],
           color='r', s=50, label="Static Points")
ax.plot(sol.y[0], sol.y[1], sol.y[2], color='b', linewidth=1.5)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title(f"a = {a}, b = {b}")

plt.show()
