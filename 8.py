from math import *
from scipy import integrate as integ
from sympy import diff
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
    for k in X:
        s *= (x - k)
    return s


I = integ.quad(F, a, b)[0]
print("Точное значение интеграла", I)
MU = np.zeros(2*N)
for i in range(2*N):
    MU[i] = integ.quad(lambda x: p(x) * x ** i, a, b)[0]
X1 = np.zeros([N, N])
for i in range(N):
    for j in range(N):
        X1[i][j] = MU[j+i]
X2 = np.zeros(N)
for i in range(N):
    X2[i] = -MU[N+i]

A = np.linalg.solve(X1, X2)
z = A.tolist()
z.append(1)
rez = np.roots(z[::-1])
# rez = [i.real for i in rez]
XZ = np.zeros([N, N])
for i in range(N):
    for j in range(N):
        XZ[i][j] = rez[j] ** i
MU_N = MU[:N]
A2 = np.linalg.solve(XZ, MU_N)
S = 0
for i in range(N):
    S += A2[i] * f(rez[i])
print("Вычисленное начение интеграла", S)
print("Разность точного интеграла и вычисленного", abs(I - S))
expr = "2 * cos(2.5 * x) * exp(x / 3) + 4 * sin(3.5 * x) * exp(-3 * x) + x"
for i in range(2*N):
    expr = diff(expr)
dfx_e = expr.evalf(subs={"x": N//2})
dd = []
for i in X:
    dd.append(abs(expr.evalf(subs={"x": i})))
M = max(dd)
R = (M / factorial(2*N)) * integ.quad(lambda x: abs(p(x) * omega(x) ** 2), a, b)[0]
print("Значение оценочной методической погрешности ", R)
