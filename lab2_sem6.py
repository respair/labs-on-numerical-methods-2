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
    if t < 1:
        return 1/k/s
    else:
        return 0



def u0n(x):
    return 0 #-x*q/(k*s)+l/(k*s)


alpha = [0, 0]
beta = [1, 1]


def koef(old, time, tau, xrange, h):
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
        d[i] = (tau * a_**2 / h ** 2)*(k-1)*(old[i+1] - 2*old[i] + old[i-1]) - tau * func(xrange[i], (time - 0.5)*tau) \
               - old[i]

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
    tn = 4
    n = 50
    
    h = (right_end - left_end) / n
    tau = h / 2
    
    time_range = np.linspace(t0, tn, int((tn - t0) / tau))
    xrange = np.arange(left_end, right_end, h)

    next = np.zeros(len(xrange))
    old_next = np.zeros(len(xrange))#u0n(xrange) #тот что был до некст

    for j in range(1, len(time_range) + 1):
        next = koef(old_next, j, tau, xrange, h)
        old_next = next

    plt.xlabel("x")
    plt.ylabel("u")
    plt.plot(xrange, next, color='r', label="t = 4")
    plt.grid()
    plt.legend()

    plt.show()


graph()
