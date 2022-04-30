import numpy as np
import matplotlib.pyplot as plt

a_ = 1
q = 1
l = 5
k = 1
s = 1


def func(x, t):
    return 0


def gamma1(t):
    return 0


def gamma0(t):
    return 0


def u0n(x):
    return -x*q/(k*s)+l/(k*s)


alpha = [1, 0]
beta = [0, 1]


def koef(prev, time, tau, xrange, h):
    d = np.zeros(len(xrange))
    d[0] = gamma0(time*tau)
    d[-1] = gamma1(time*tau)

    b = np.zeros(len(xrange))
    b[0] = - alpha[0]/h + beta[0]
    b[-1] = alpha[1]/h + beta[1]

    c = np.zeros(len(xrange))
    c[0] = alpha[0]/h

    a = np.zeros(len(xrange))
    a[-1] = -alpha[1]/h

    for i in range(1, len(xrange)-1):
        a[i] = 1/2 * (tau * a_**2 / h**2)
        b[i] = -1 - 2 * 1/2 * (tau * a_**2 / h**2)
        c[i] = 1/2 * (tau * a_**2 / h**2)
        d[i] = (tau * a_**2 / h ** 2)*(k-1)*(prev[i+1] - 2*prev[i] + prev[i-1]) - tau * func(xrange[i], (time - 0.5)*tau) \
               - prev[i]

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
    t0 = 1
    tn = 5
    n = 50
    h = (right_end - left_end) / n
    tau = h / 2
    time_range = np.linspace(t0, tn, int((tn - t0) / tau))
    xrange = np.arange(left_end, right_end, h)

    prev = u0n(xrange)
    next = np.zeros(len(xrange))

    for i in range(1, len(time_range) + 1):
        next = koef(prev, i, tau, xrange, h)
        prev = next

    plt.xlabel("x")
    plt.ylabel("u")
    plt.grid()
    plt.plot(xrange, next, color='r', label="Численное решение")
    plt.legend()

    plt.show()


graph()
