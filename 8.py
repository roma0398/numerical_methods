from math import factorial
from sympy import (
    integrate,
    symbols,
    diff,
    cos,
    exp,
    sin,
    simplify,
    expand,
    Function,
    S,
    Rational,
    solve,sympify

)
import numpy as np

a = Rational('1.5')
b = Rational('3.3')
a_ = 1.5
b_ = 3.3
alpha = S(1) / 3
beta = 0
N = 3

x = symbols("x")
f = 2 * cos(2.5 * x) * exp(x * alpha) + 4 * sin(3.5 * x) * exp(-3 * x) + x
p = (x - a) ** -alpha * (b - x) ** -beta
fun = f * p
MU = np.zeros(2*N)
for i in range(2*N):
    MU[i] = float(integrate(p*x**i, (x, a, b)))
X1 = np.zeros([N, N])
for i in range(N):
    for j in range(N):
        X1[i][j] = MU[j+i]
X2 = np.zeros(N)
for i in range(0, N):
    X2[i] = -MU[N+i]
M = np.linalg.solve(X1, X2)
omega = x ** N
for i in np.arange(N-1, -1, -1):
    omega += M[i]*x**i
rez = []
xx = solve(omega)
print(type(xx[0]))
def f(k):
    return k.evalf().real
for i in xx:
    rez.append(f(i))
X = np.zeros([N, N])
for i in range(N):
    for j in range(N):
        X[i][j] = rez[j] ** i
MU_N = MU[:N]
A = np.linalg.solve(X, MU_N)
S = 0
for i in range(N):
    S = S+A[i]*f.evalf(subs={'x': rez[i]})
I = integrate(fun, (x, a, b))
print(S)
print(I)
print(S-I)
dfn = diff(f, (x, 2*N))
x_ar = np.arange(a, b, 0.001)
Mn = max(abs(dfn.evalf(subs={'x': x_ar})))
omega2 = 1
for i in range(N):
    omega2 *= (x - rez(i))
f1 = p * omega2 * omega2
f2 = abs(f1)
R = Mn/factorial(2*N)*integrate(f2, (x, a, b))
print(R)