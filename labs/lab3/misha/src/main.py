from lab3 import lib
import numpy as np
import cffi
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

ffi = cffi.FFI()

# Инициализация данных
N = 1000
count = int(1e6)
tau = 0.01
alpha = 0
beta = 305

q = np.zeros(N, dtype=np.float64)
v = np.zeros(N, dtype=np.float64)
a = np.zeros(N, dtype=np.float64)
q[N // 2] = -0.5
q[N // 2 - 1] = 0.5

q_ptr = ffi.cast("double*", ffi.from_buffer(q))
v_ptr = ffi.cast("double*", ffi.from_buffer(v))
a_ptr = ffi.cast("double*", ffi.from_buffer(a))

# lib.SimplexVerle(q_ptr, v_ptr, N, 0.0, 100.0, 0.01, 1000000, 1.0)

data = []
print(f'Начальный гамильтон: {lib.H(q_ptr, v_ptr, N, 1, alpha, beta)}')
for i in range(count):
    lib.SimplexVerleNew(q_ptr, v_ptr, a_ptr, N, alpha, beta, tau, 1)
    if i % 1000 == 0:
        data.append(v.copy())
print(f'Конечный гамильтон: {lib.H(q_ptr, v_ptr, N, 1, alpha, beta)}')

fig = plt.figure()
ax = plt.axes(xlim=(0, len(data[0])), ylim=(-1, 1))

x = np.arange(len(data[0]))
line, = ax.plot(x, data[0], lw=2)

st = 0
fr = len(data)


def init():
    line.set_data(x, data[0])
    return line,


def animate(i):
    if st + i < fr:
        line.set_data(x, data[st + i])
        return line,
    sys.exit(1)


anim = FuncAnimation(fig, animate, frames=1000, interval=50, blit=True)
# anim.save('2_fpu_betta.mp4', writer='ffmpeg', fps=25)

plt.show()

# 165 0.5 один
#
#
