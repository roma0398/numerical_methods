from math import *
import scipy.integrate as spint
from sympy import diff
import numpy as np

a = 1.5
b = 3.3
alpha = 1 / 3
beta = 0


def f(x):
    return 2 * cos(2.5 * x) * exp(x * alpha) + 4 * sin(3.5 * x) * exp(-3 * x) + x


def p(x):
    return 1 / ((x - a) ** alpha * (b - x) ** beta)


def F(x):
    return f(x) * p(x)


def px(x):
    return x / ((x - a) ** alpha * (b - x) ** beta)


def pxx(x):
    return x ** 2 / ((x - a) ** alpha * (b - x) ** beta)


def omega(x):
    return (x - x1) * (x - x2) * (x - x3)


def dfx(x):
    return expr.evalf(subs={"x": x})


def p_omega_dfx(x):
    return p(x) * omega(x) * dfx_e


def abs_p_omega(x):
    return abs(p(x) * omega(x))


x1 = a
x2 = (a + b) / 2
x3 = b
I = spint.quad(F, a, b)[0]
print("Точное значение интеграла", I)

mu1 = spint.quad(p, a, b)[0]
mu2 = spint.quad(px, a, b)[0]
mu3 = spint.quad(pxx, a, b)[0]
M1 = np.array([[1, 1, 1], [x1, x2, x3], [x1 ** 2, x2 ** 2, x3 ** 2]])
v1 = np.array([mu1, mu2, mu3])
A = np.linalg.solve(M1, v1)
X = [x1, x2, x3]
sum = 0
for i in range(3):
    sum += A[i] * f(X[i])
print("Вычисленное начение интеграла", sum)
expr = diff(
    diff(diff("2 * cos(2.5 * x) * exp(x / 3) + 4 * sin(3.5 * x) * exp(-3 * x) + x"))
)
dfx_e = dfx(x2)
R_1 = spint.quad(p_omega_dfx, a, b)[0] / factorial(3)
print("Значение точной методической погрешности", R_1)
x_area = np.linspace(a, b, 2000)
dd = []
for i in x_area:
    dd.append(abs(dfx(i)))
M = max(dd)
R = (M / factorial(3)) * spint.quad(abs_p_omega, a, b)[0]
print("Значение оценочной методической погрешности ", R)
print("Разность точного интеграла и вычисленного", I - sum)
