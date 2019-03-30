from math import *
from scipy import integrate as integ
from sympy import diff
import numpy as np


def NK(z1, z2, N):
    X = np.linspace(z1, z2, N)
    MU = []
    for i in range(len(X)):
        MU.append(integ.quad(lambda x: p(x) * x ** i, z1, z2)[0])
    M1 = []
    for i in range(len(X)):
        M1.append([j ** i for j in X])
    v1 = MU
    A = np.linalg.solve(M1, v1)
    ssum = 0
    for i in range(N):
        ssum += A[i] * f(X[i])
    return ssum


a = 1.5
b = 3.3
alpha = 1 / 3
beta = 0
k = 50
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
    S += NK(z_1, z_2, N)

print("\nКвадратурная формула Ньютона-Котеса")
print("Вычисленное начение интеграла", S)
print("Длина шага", h)
print("Разность точного интеграла и вычисленного", abs(I - S))
