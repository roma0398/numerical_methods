from math import *
from scipy import integrate as integ
from sympy import diff
import numpy as np

a = 1.5
b = 3.3
alpha = 1 / 3
beta = 0
k = 4
h = (b - a) / k
K = np.linspace(a, b, k)


def f(x):
    return 2 * cos(2.5 * x) * exp(x * alpha) + 4 * sin(3.5 * x) * exp(-3 * x) + x


def p(x):
    return 1 / ((x - a) ** alpha * (b - x) ** beta)


def F(x):
    return f(x) * p(x)


def omega(x):
    s = 1
    for w in X:
        s *= (x - w)
    return s


I = integ.quad(F, a, b)[0]
print("Точное значение интеграла", I)
N = 3
S = 0
for i in range(k):
    z_1 = a + i * h
    z_2 = a + (i+1)*h
    X = np.linspace(z_1, z_2, N)
    MU = np.zeros(2 * N)
    for i in range(2 * N):
        MU[i] = integ.quad(lambda x: p(x) * x ** i, z_1, z_2)[0]
    X1 = np.zeros([N, N])
    for i in range(N):
        for j in range(N):
            X1[i][j] = MU[j + i]
    X2 = np.zeros(N)
    for i in range(N):
        X2[i] = -MU[N + i]
    A = np.linalg.solve(X1, X2)
    z = A.tolist()
    z.append(1)
    rez = np.roots(z[::-1])
    rez = [i.real for i in rez]
    XZ = np.zeros([N, N])
    for i in range(N):
        for j in range(N):
            XZ[i][j] = rez[j] ** i
    MU_N = MU[:N]
    A2 = np.linalg.solve(XZ, MU_N)
    ssum = 0
    for i in range(N):
        ssum += A2[i] * f(rez[i])
    S += ssum

print("\nКвадратурная формула Гаусса")
print("Вычисленное начение интеграла", S)
print("Длина шага", h)
print("Разность точного интеграла и вычисленного", abs(I - S))