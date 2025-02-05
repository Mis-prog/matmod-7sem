from lab3 import lib
import numpy as np
import cffi
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
from matplotlib.widgets import Button, Slider

ffi = cffi.FFI()

# Инициализация данных
N = 1000
count = 1000000
tau = 0.01
alpha = 0
beta = 200
q = np.zeros(N, dtype=np.float64)
v = np.zeros(N, dtype=np.float64)
a = np.zeros(N, dtype=np.float64)
q[N // 2] = 0.5
q[N // 2 - 1] = -0.5

q_ptr = ffi.cast("double*", ffi.from_buffer(q))
v_ptr = ffi.cast("double*", ffi.from_buffer(v))
a_ptr = ffi.cast("double*", ffi.from_buffer(a))


# lib.SimplexVerle(q_ptr, v_ptr, N, 0.0, 100.0, 0.01, 1000000, 1.0)

data = []
for i in range(count):
    lib.SimplexVerleNew(q_ptr, v_ptr, a_ptr, N, alpha, beta, tau, 1)
    if i % 1000 == 0:
        data.append(v.copy())

fig = plt.figure()
ax = plt.axes(xlim=(0, len(data[0])), ylim=(-4, 4))

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


anim = FuncAnimation(fig, animate, frames=fr, interval=10, blit=True)

# --- Кнопки ---
def toggle_animation(event):
    global anim_running
    if anim_running:
        anim.event_source.stop()
    else:
        anim.event_source.start()
    anim_running = not anim_running


def reset_simulation(event):
    global st
    st = 0
    line.set_data(x, data[0])
    fig.canvas.draw()


def update_beta(val):
    """Обновление beta"""
    global beta
    beta = beta_slider.val


# --- Размещение кнопок ---
ax_start = plt.axes([0.7, 0.05, 0.1, 0.075])
ax_reset = plt.axes([0.81, 0.05, 0.1, 0.075])
ax_beta = plt.axes([0.25, 0.05, 0.4, 0.03])

btn_start = Button(ax_start, "Старт/Стоп")
btn_reset = Button(ax_reset, "Сброс")
beta_slider = Slider(ax_beta, "β", 50, 500, valinit=beta)

btn_start.on_clicked(toggle_animation)
btn_reset.on_clicked(reset_simulation)
beta_slider.on_changed(update_beta)

plt.show()