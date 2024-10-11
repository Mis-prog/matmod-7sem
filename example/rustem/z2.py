from scipy import integrate
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib import animation
import random
from matplotlib.animation import FuncAnimation
import math
import numpy as np

m1, m2, m3 = 2 * (10 ** 30), 1.9 * (10 ** 27), 1.5 * (10 ** 23)
G = 6.67 * (10 ** -11)
m0 = 75
mt = 4950
T = 1600
koef = 0.001
u = 3170
Rpl = 71500 * (10 ** 3)
H = 900 * (10 ** 3)
Rsat = 2634 * (10 ** 3)


def m(t):
    if t >= T:
        return m0
    else:
        return (m0 + mt) / (1 - koef) - mt * t / T


def dm(t):
    if t >= T:
        return 0
    else:
        return -mt / T


r12x = 0
r12y = 0


def event2(t, y):
    if not hasattr(event2, "counter"):
        event2.stop = -1
    if event2.stop == -1 and (math.sqrt((y[2] - y[0]) ** 2 + (y[3] - y[1]) ** 2) - Rsat) < 0:
        event2.stop = t
    return t - event2.stop


event2.terminal = True


def F2(t, y):
    rx, ry, r13x, r13y, vx, vy, v3x, v3y = y
    r = math.sqrt(rx ** 2 + ry ** 2)
    v = math.sqrt(vx ** 2 + vy ** 2)
    r2 = math.sqrt((rx - r12x) ** 2 + (ry - r12y) ** 2)
    r3 = math.sqrt((rx - r13x) ** 2 + (ry - r13y) ** 2)
    r13 = math.sqrt(r13x ** 2 + r13y ** 2)
    r23x = (r13x - r12x)
    r23y = (r13y - r12y)
    r23 = math.sqrt((r13x - r12x) ** 2 + (r13y - r12y) ** 2)

    return [vx, vy, v3x, v3y,
            (-u * dm(t) / m(t) * vx / v + G * (
                    -m1 * rx / (r ** 3) - m2 * (rx - r12x) / (r2 ** 3) - m3 * (rx - r13x) / (r3 ** 3))),
            (-u * dm(t) / m(t) * vy / v + G * (
                    -m1 * ry / (r ** 3) - m2 * (ry - r12y) / (r2 ** 3) - m3 * (ry - r13y) / (r3 ** 3))),
            G * (-(m1 / (r13 ** 3) * r13x) - (m2 / (r23 ** 3) * r23x)),
            G * (-(m1 / (r13 ** 3) * r13y) - (m2 / (r23 ** 3) * r23y))
            ]


def z2(r13x0, r13y0, v3x0, v3y0, tau, r12x0, r12y0, fi, m):
    global r12x, r12y
    global mt
    print(m)
    mt = m
    r12x = r12x0
    r12y = r12y0
    r3x = r13x0 - r12x
    r3y = r13y0 - r12y
    r3 = math.sqrt(r3x ** 2 + r3y ** 2)
    v0 = 1.0 * math.sqrt(G * m2 / (Rpl + H))
    rx0 = (Rpl + H) * (r3x * math.cos(fi) - r3y * math.sin(fi)) / r3
    ry0 = (Rpl + H) * (r3x * math.sin(fi) + r3y * math.cos(fi)) / r3
    r0 = math.sqrt(rx0 ** 2 + ry0 ** 2)
    vx0 = -v0 * ry0 / r0
    vy0 = v0 * rx0 / r0
    rx0 = r12x + rx0
    ry0 = r12y + ry0
    s0 = [rx0, ry0, r13x0, r13y0, vx0, vy0, v3x0, v3y0]
    t_span = (0, 800000)
    satel = solve_ivp(F2, t_span, s0, max_step=tau, atol=1, rtol=1, method="RK45", events=event2)
    return satel


def optimization(solution, m):
    R = [math.sqrt((solution.y[2][i] - solution.y[0][i]) ** 2 + (solution.y[3][i] - solution.y[1][i]) ** 2) - Rsat for i
         in
         range(len(solution.y[0]))]
    Res = math.sqrt((solution.y[2][-1] - solution.y[0][-1]) ** 2 + (solution.y[3][-1] - solution.y[1][-1]) ** 2) - Rsat
    if Res < 0:
        Res = 0
    Res += m
    return max(min(R), 0) + m


def show_animation(satel, r13x0, r13y0, tau):
    r13 = math.sqrt(r13x0 ** 2 + r13y0 ** 2)
    alpha = math.acos(r13x0 / r13)
    s = int(300 / tau)
    kr = 10000000
    sx = [(satel.y[0][s * i] - r12x) * math.cos(alpha) - (satel.y[1][s * i] - r12y) * math.sin(alpha) for i in
          range(int(len(satel.y[0]) / s))]
    sy = [(satel.y[0][s * i] - r12x) * math.sin(alpha) + (satel.y[1][s * i] - r12y) * math.cos(alpha) for i in
          range(int(len(satel.y[0]) / s))]
    lx = [(satel.y[2][s * i] - r12x) * math.cos(alpha) - (satel.y[3][s * i] - r12y) * math.sin(alpha) for i in
          range(int(len(satel.y[0]) / s))]
    ly = [(satel.y[2][s * i] - r12x) * math.sin(alpha) + (satel.y[3][s * i] - r12y) * math.cos(alpha) for i in
          range(int(len(satel.y[0]) / s))]
    fig, ax = plt.subplots()
    line, = ax.plot(sx, sy, color='g', label='ракета')
    line2, = ax.plot(lx, ly, label='спутник')
    plt.plot(0, 0, 'ro')
    plt.legend()
    moon = plt.Circle((r13x0 - r12x, r13y0 - r12y), Rsat, color='r')
    rocket = plt.Circle((r13x0 - r12x, r13y0 - r12y), 5 * Rsat, color='c')
    p = ax.add_patch(moon)
    q = ax.add_patch(rocket)

    def animate(i):
        line.set_xdata(sx[0:i])
        line.set_ydata(sy[0:i])  # update the data
        line2.set_xdata(lx[0:i])  # update the data
        line2.set_ydata(ly[0:i])  # update the data
        moon.set_center((lx[i], ly[i]))
        rocket.set_center((sx[i], sy[i]))
        return moon, rocket, line, line2,

    ani = animation.FuncAnimation(fig, animate, np.arange(1, int(len(satel.y[0]) / s)), interval=1)
    plt.show()
    ani.save('myAnimation.gif', writer='pillow', fps=30)
    print('gif saved')


def NelderMid():
    a1 = random.random() * math.pi * 2
    mt1 = random.randint(1, 20000)
    tau = 0.1
    a2 = random.random() * math.pi * 2
    mt2 = random.randint(1, 20000)
    a3 = random.random() * math.pi * 2
    mt3 = random.randint(1, 20000)
    res1 = optimization(z2(-166864782814.52066, -751687780697.0953, 10608.237381615552, -2321.8734688953678, tau,
                           -166632056766.0584, -750639289651.4109, a1, mt1), mt1)
    res2 = optimization(z2(-166864782814.52066, -751687780697.0953, 10608.237381615552, -2321.8734688953678, tau,
                           -166632056766.0584, -750639289651.4109, a2, mt2), mt2)
    res3 = optimization(z2(-166864782814.52066, -751687780697.0953, 10608.237381615552, -2321.8734688953678, tau,
                           -166632056766.0584, -750639289651.4109, a3, mt3), mt3)

    for i in range(15):
        if res2 == max([res1, res2, res3]):
            res2, res1 = res1, res2
            a2, a1 = a1, a2
            mt2, mt1 = mt1, mt2
            if res2 < res3:
                res2, res3 = res3, res2
                a2, a3 = a3, a2
                mt2, mt3 = mt3, mt2
        if res3 == max([res1, res2, res3]):
            res3, res1 = res1, res3
            a3, a1 = a1, a3
            mt3, mt1 = mt1, mt3
            if res2 < res3:
                res2, res3 = res3, res2
                a2, a3 = a3, a2
                mt2, mt3 = mt3, mt2
        a0 = (a2 + a3) / 2
        mt0 = (mt2 + mt3) / 2
        an = (a2 + a3) - a1
        mtn = (mt2 + mt3) - mt1
        an = max(0, an)
        mtn = max(0, mtn)
        res = optimization(z2(-166864782814.52066, -751687780697.0953, 10608.237381615552, -2321.8734688953678, tau,
                              -166632056766.0584, -750639289651.4109, an, mtn), mtn)
        if res < max(res2, res3):
            an2 = a0 + 2 * (an - a0)
            mtn2 = mt0 + 2 * (mtn - mt0)
            an2, mtn2 = max(0, an), max(0, mtn)
            res4 = optimization(
                z2(-166864782814.52066, -751687780697.0953, 10608.237381615552, -2321.8734688953678, tau,
                   -166632056766.0584, -750639289651.4109, an2, mtn2), mtn2)
            if res4 < min(res2, res3):
                an, mtn, res = an2, mtn2, res4
        if res > max(res3, res2):
            an = a0 + 0.5 * (a1 - a0)
            mtn = mt0 + 0.5 * (mt1 - mt0)
            an, mtn = max(0, an), max(0, mtn)
            res = z2(-166864782814.52066, -751687780697.0953, 10608.237381615552, -2321.8734688953678, tau,
                     -166632056766.0584, -750639289651.4109, an, mtn)
        if res < max(res3, res2):
            a1, mt1, res1 = an, mtn, res
        else:
            a1 = a3 + 0.5 * (a1 - a3)
            mt1 = mt3 + 0.5 * (mt1 - mt3)
            a2 = a3 + 0.5 * (a2 - a3)
            mt2 = mt3 + 0.5 * (mt2 - mt3)
            a1 = max(0, a1)
            mt1 = max(0, mt1)
            a2 = max(0, a2)
            mt2 = max(0, mt2)
            res1 = z2(-166864782814.52066, -751687780697.0953, 10608.237381615552, -2321.8734688953678, tau,
                      -166632056766.0584, -750639289651.4109, a1, mt1)
            res2 = z2(-166864782814.52066, -751687780697.0953, 10608.237381615552, -2321.8734688953678, tau,
                      -166632056766.0584, -750639289651.4109, a2, mt2)
