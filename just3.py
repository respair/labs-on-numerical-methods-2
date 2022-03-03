import numpy as np
import math
import matplotlib.pyplot as plt


def u_0(x, t):
    return x ** 2 + np.arccos(x * t / 2)


def u_x0(x):
    return x ** 2 + math.pi / 2


def ut_x0(x):
    return -x / 2


def func(x, t):
    return -1 + x * 0.5 * t * (t ** 2 - 2 * (x ** 2)) / (4 - (x ** 2) * t ** 2) ** (3 / 2)


# граничные условия
# a[0] * u(0, t) + b[0] * ux(0, t) = g[0]
# a[1] * u(1, t) + b[1] * ux(1, t) = g[1]
alpha = [0, 1]
beta = [1, 1]


def gamma0(t):
    return -t / 2


def gamma1(t):
    return 3 + np.arccos(t / 2) - t / (4 - t ** 2) ** (1 / 2)


def unext_2(u_, a, h, t, now):
    u = np.zeros(len(u_[0]))
    for i in range(1, len(u) - 1):
        u[i] = 2 * u_[1, i] - u_[0, i] + ((a * t / h) ** 2) * \
               (u_[1, i + 1] - 2 * u_[1, i] + u_[1, i - 1]) + t ** 2 * func(i * h, now - t)
    u[0] = (gamma0(now) + ((beta[0] / (2 * h)) * (u[2] - 4 * u[1]))) / (alpha[0] - (3 * beta[0]) / (2 * h))
    u[-1] = (gamma1(now) + (beta[1] / (2 * h)) * (-4 * u[-2] + u[-3])) / (alpha[1] - (3 * beta[1]) / (2 * h))
    return u


def u0n_2(u, t, arr_x, a):
    return u + t * ut_x0(arr_x) + ((t ** 2) / 2) * (a ** 2 * u_x0(arr_x) + func(arr_x, 0))


def graph():
    C = 1 / 2
    a = np.sqrt(1 / 2)

    x0 = 0
    xn = 1

    h0 = 0.001
    hn = 0.01
    step = 0.001
    h = np.arange(h0, hn, step)

    t0 = 0
    T = 20

    error = np.zeros(len(h))
    i = 0

    while h0 <= hn:
        N = round((xn - x0) / h0)
        arr_x = np.linspace(x0, xn, N)
        t = C * h0 / a
        tn = t0 + t * T
        arr_t = np.linspace(t0 + 3 * t, tn, T)
        u = np.zeros((3, N))
        u[0] = u_x0(arr_x)
        u[1] = u0n_2(u[0], t, arr_x, a)
        u[2] = unext_2(u[0:2], a, h0, t, t0 + 2 * t)
        for t_ in arr_t:
            u[0] = u[1]
            u[1] = u[2]
            u[2] = unext_2(u[0:2], a, h0, t, t_)
        error[i] = np.max(np.abs(u[2] - u_0(arr_x, tn)))
        h0 += step
        i += 1

    h = np.log(h)
    error = np.log(error)

    plt.suptitle('Зависимость логарифма абсолютной погрешности от логарифма шага интегрирования')
    plt.subplot(1, 1, 1)
    plt.xlabel("log(h)")
    plt.ylabel("log(max(|Δu|))")
    plt.grid()
    plt.plot(h, error, color='k')

    plt.show()


graph()
