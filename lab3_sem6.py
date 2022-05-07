import numpy as np
import matplotlib.pyplot as plt

#a_ = 1
q = 1
l = 2
k = 1
s = 1
n = 1
u0_ = 1


def lambda_(u):
    return u ** n


def alpha_():
    return 1 / (n ** (1 / 2))


def func(x, t, lend):
    if x == lend:
        return 0
    if x == 0 and t != 0:
        return 1
    if t == 0:
        return 0
    else:
        z = x / (u0_ ** n * t ** (k * n + 1)) ** 0.5
        if z <= 1 / (n ** 2):
            return alpha_() * n * (alpha_() - z)


def gamma1(t):
    return 0


def gamma0(t):
    return u0_ * t ** k


def u0n(x):
    return 0  # -x*q/(k*s)+l/(k*s)


alpha = [0, 0]
beta = [1, 1]


def u0(x, t, xn):
    if x <= (alpha_() * (u0_ ** n * t ** (k * n + 1)) ** 0.5):
        return u0_ * t ** k * (func(x, t, xn) ** (1 / n))
    else:
        return 0


def koef(old, u, time, tau, xrange, h):
    d = np.zeros(len(xrange))
    d[0] = gamma0(time * tau)
    d[-1] = gamma1(time * tau)

    b = np.zeros(len(xrange))
    b[0] = - alpha[0] / h + beta[0]
    b[-1] = alpha[1] / h + beta[1]

    c = np.zeros(len(xrange))
    c[0] = alpha[0] / h

    a = np.zeros(len(xrange))
    a[-1] = -alpha[1] / h

    for i in range(1, len(xrange) - 1):
        a[i] = tau / (2 * h ** 2) * (lambda_(u[i]) + lambda_(u[i - 1]))
        b[i] = -tau / (2 * h ** 2) * (
                    2 * lambda_(u[i]) + lambda_(u[i - 1]) + lambda_(u[i + 1])) - 1
        c[i] = tau / (2 * h ** 2) * (lambda_(u[i]) + lambda_(u[i + 1]))
        d[i] = -old[i]

    return solve(a, b, c, d, xrange)


def solve(a, b, c, d, x):
    A = np.zeros(len(x))
    A[0] = - c[0] / b[0]

    B = np.zeros(len(x))
    B[0] = d[0] / b[0]

    for i in range(1, len(x)):
        A[i] = -c[i] / (b[i] + a[i] * A[i - 1])
        B[i] = (d[i] - a[i] * B[i - 1]) / (b[i] + a[i] * A[i - 1])

    result = np.zeros(len(x))
    result[-1] = B[-1]

    for i in range(len(x) - 2, -1, -1):
        result[i] = B[i] + A[i] * result[i + 1]

    return result


def graph():
    left_end = 0
    right_end = l
    t0 = 0
    tn = 2
    n = 60
    h = (right_end - left_end) / n
    tau = 0.01
    time_range = np.linspace(t0, tn, int((tn - t0) / tau))
    xrange = np.arange(left_end, right_end, h)

    old_next = np.zeros(len(xrange))  # u0n(xrange) #тот что был до некст
    u_ = np.zeros(len(xrange))

    for i in range(1, len(time_range)):
        for j in range(5):
            u_ = koef(old_next, u_, i, tau, xrange, h)
        old_next = u_
    next = u_

    u = np.zeros(len(xrange))
    for i in range(len(xrange)):
        el = xrange[i]
        u[i] = u0(el, tn, xrange[-1])

    err = abs(u - next)

    plt.figure("Решение")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.plot(xrange, next, label="численное решение")
    plt.plot(xrange, u, color='g', label="аналитическое решение")
    plt.grid()
    plt.legend()

    plt.figure("Ошибка")
    plt.xlabel("x")
    plt.ylabel("|delta|")
    plt.plot(xrange, err, label="ошибка")
    plt.grid()
    plt.legend()

    plt.show()


graph()
