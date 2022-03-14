import numpy as np
import math
import pylab


def u_0(x, t):
    return x ** 2 + np.arccos(x * t / 2)


def u_x0(x):
    return x ** 2 + math.pi / 2


def ut_x0(x):
    return -x / 2


def func(x, t):
    return -1 + x * 0.5 * t * (t ** 2 - 2 * (x ** 2)) / (4 - (x ** 2) * t ** 2) ** (3 / 2)


# граничные условия
# a[0] * u(0, t) + b[0] * ux(0, t) = c[0]
# a[1] * u(1, t) + b[1] * ux(1, t) = c[1]
a_ = [0, 1]
b = [1, 1]


def c0(t):
    return -t / 2


def c1(t):
    return 3 + np.arccos(t / 2) - t / (4 - t ** 2) ** (1 / 2)


def unext(indicator, u_, a, h, t, now):
    l = len(u_[0])
    u = np.zeros(l)
    i = 1
    if indicator == 2:
        while i < l - 1:
            u[i] = 2 * u_[1][i] - u_[0][i] + ((a * t / h) ** 2) * \
                   (u_[1][i + 1] - 2 * u_[1][i] + u_[1][i - 1]) + t ** 2 * func(i * h, now - t)
            i += 1
        u[0] = (c0(now) + ((b[0] / (2 * h)) * (u[2] - 4 * u[1]))) / (a_[0] - (3 * b[0]) / (2 * h))
        u[-1] = (c1(now) + (b[1] / (2 * h)) * (-4 * u[-2] + u[-3])) / (a_[1] - (3 * b[1]) / (2 * h))
        return u
    else:
        while i < l - 1:
            u[i] = 2 * u_[1][i] - u_[0][i] + ((a * t / h) ** 2) * \
                   (u_[1][i + 1] - 2 * u_[1][i] + u_[1][i - 1]) + t ** 2 * func(i * h, now - t)
            i += 1
        u[0] = (c0(now) - b[0] * u[1] / h) / (a_[0] - b[0] / h)
        u[-1] = (c1(now) + b[1] * u[-2] / h) / (a_[1] + b[1] / h)
        return u


def u0n(indicator, u, t, arr_x, a):
    if indicator == 2:
        return u + t * ut_x0(arr_x) + ((t ** 2) / 2) * (a ** 2 * u_x0(arr_x) + func(arr_x, 0))
    else:
        return u + t * ut_x0(arr_x)


def graph():
    indicator = 2
    C = 1 / 3
    a = np.sqrt(1 / 2)

    x0 = 0
    xn = 1

    h0 = 0.001
    hn = 0.01
    step = 0.001
    h = np.arange(h0, hn, step)
    hh = []

    t0 = 0
    # tn = 0.5
    T = 20 # количество точек

    err = np.zeros(len(h))
    i = 0

    for h0 in h:
        N = round((xn - x0) / h0)
        u = [np.zeros(N), np.zeros(N), np.zeros(N)]
        arr_x = np.linspace(x0, xn, N)  # np.arange(x0,xn,step)
        t = C * h0 / a
     #   T = round((tn - t0) / t)
        tn = t0 + t * T
     #   print(t)
     #   print(T)
        arr_t = np.linspace(t0 + 3 * t, tn, T)
        u[0] = u_x0(arr_x)
        u[1] = u0n(indicator, u[0], t, arr_x, a)
        temp = t0 + 2*t
        u[2] = unext(indicator, u, a, h0, t, temp)
        for t_ in arr_t:
            u[0] = u[1]
            u[1] = u[2]
            u[2] = unext(indicator, u, a, h0, t, t_)
        err[i] = np.max(np.abs(u[2] - u_0(arr_x, tn)))
     #   hh.append(h0)
     #   h0 += step
        i += 1

    pylab.plot(np.log(h), np.log(err), label=str(indicator) + " порядок")
    pylab.grid()
    pylab.legend()
    pylab.show()


graph()
