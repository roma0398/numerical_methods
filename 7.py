from math import *
from scipy import integrate as integ
from sympy import *
import numpy as np

a = 1.5
b = 3.3
alpha = 1 / 3
beta = 0
N = 5
X = np.linspace(a, b, N)


def f(x):
    return 2 * cos(2.5 * x) * exp(x * alpha) + 4 * sin(3.5 * x) * exp(-3 * x) + x


def p(x):
    return 1 / ((x - a) ** alpha * (b - x) ** beta)


def F(x):
    return f(x) * p(x)


def omega(x):
    s = 1
    for i in X:
        s *= (x - i)
    return s


I = integ.quad(F, a, b)[0]
print("Точное значение интеграла", I)


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


print("Вычисленное начение интеграла", NK(a, b, N))
print("Разность точного интеграла и вычисленного", abs(I - NK(a, b, N)))
expr = "2 * cos(2.5 * x) * exp(x / 3) + 4 * sin(3.5 * x) * exp(-3 * x) + x"
for i in range(N):
    expr = diff(expr)
dfx_e = expr.evalf(subs={"x": N//2})
dd = []
for i in X:
    dd.append(abs(expr.evalf(subs={"x": i})))
M = max(dd)
R = (M / factorial(N)) * integ.quad(lambda x: abs(p(x) * omega(x)), a, b)[0]
print("Значение оценочной методической погрешности ", R)
