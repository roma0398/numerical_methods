from math import factorial
from sympy import integrate, symbols, diff, cos, exp, sin, simplify, expand, Function, S, Rational
import numpy as np

a = Rational('1.5')
b = Rational('3.3')
a_ = 1.5
b_ = 3.3
alpha = S(1) / 3
beta = 0
N = 3
X = np.linspace(a_, b_, N)

x = symbols("x")
f = 2 * cos(2.5 * x) * exp(x * alpha) + 4 * sin(3.5 * x) * exp(-3 * x) + x
p = 1 / ((x - a) ** alpha * (b - x) ** beta)
F = f * p
d = p * x
print(integrate(x*p, (x, a, b)))
MU = np.zeros(N)
for i in range(N - 1):
    MU[i] = integrate(p * x ** i, (x, a, b))

XX = np.eye(N, N)
for i in range(N):
    for j in range(N):
        XX[i][j] = X[j] ** (i - 1)

A = np.linalg.solve(XX, MU)
sum = 0
for i in range(N):
    sum += A[i] * f.evalf(subs={'x': X[i]})
print("Вычисленное значение интеграла", sum)
dx3 = diff(f, x, 3)
omega = 1
for i in range(N):
    omega = omega * (x - X[i])

R_1 = integrate(
    p * omega * (dx3.evalf(subs={"x": X[N // 2]})), (x, a, b)) / factorial(N)


dd = []
for i in X:
    dd.append(abs(dx3.evalf(subs={"x": i})))
M = max(dd)
R = (M / factorial(N)) * integrate(abs(p * omega), (x, a, b))
print("Значение оценочной методической погрешности ", R)
I = integrate(p * f, (x, a, b))
print("Разность точного интеграла и вычисленного", I - sum)
