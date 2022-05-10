import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def a(indicator, x):
    if indicator == 2:
        return x
    else:
        return x**(1/2)


def func(indicator, x, a_):
    if indicator == 1:
        return a_*x
    if indicator == 2:
        return a_ * (x**2) / 2
    else:
        return 2/3*x**(3/2)


def gamma1(x, t):
    return 0


def gamma0(x, t):
    return 0


def u0(indicator, x):
    if indicator == 1:
        return np.heaviside(x - 1, 1) * np.heaviside(2 - x, 1)
    if indicator == 2:
        return (np.sin(np.pi * x / 2) ** 2) * np.heaviside(2 - x, 1) * np.heaviside(x, 1)
    else:
        return np.heaviside(x - 0.1, 1) * np.abs(np.sin(4*x)) * np.heaviside(1.5 - x, 1)


def Lax_Friedrichs(indicator, prev, xrange, h, t, a_):
    u = np.zeros(len(xrange))
    for i in range(1, len(xrange) - 1):
        u[i] = 1/2 * (prev[i + 1] + prev[i - 1]) - t / (2 * h) * (func(indicator, prev[i + 1], a_)
                                                                    - func(indicator, prev[i - 1], a_))
    u[0] = gamma0(xrange, t)
    u[-1] = gamma1(xrange, t)
    return u


def anim():
    indicator = 3
    a_ = 1
    start = 0
    end = 10
    t0 = 0
    tn = 4
    C = 0.98
    m = 500
    h = (end - start) / m
    xrange = np.linspace(start, end, m)

    u = np.zeros((2, m))
    u[0] = u0(indicator, xrange)

    if indicator == 1:
        t = C * h / a_
    else:
        t = (C * h) / max(a(indicator, u[0]))

    u[1] = Lax_Friedrichs(indicator, u[0], xrange, h, t, a_)

    fig = plt.figure()
    ax = plt.axes(xlim=(start, end), ylim=(-1.5, 1.5))
    line, = ax.plot([], [], lw=3)

    def init():
        line.set_data([], [])
        return line,

    def animate(i):
        if i in [0, 1]:
            u_rez = u[i]
        else:
            u[0] = u[1]
            u[1] = Lax_Friedrichs(indicator, u[0], xrange, h, t, a_)
            u_rez = u[1]

        line.set_data(xrange, u_rez)
        return line,

    anim = FuncAnimation(fig, animate, init_func=init,
                         frames=int((tn - t0) // t), interval=20, blit=True)

    anim.save('ex_test.gif', writer='imagemagick')


anim()
