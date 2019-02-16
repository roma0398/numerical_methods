from math import *
import scipy.integrate as spint

a = 1.5
b = 3.3
alpha = 1 / 3
beta = 0


def f(x):
    return 2 * cos(2.5 * x) * exp(x * alpha) + 4 * sin(3.5 * x) * exp(-3 * x) + x


def p(x):
    return (x - a) ** alpha * (b - x) ** beta


def F(x):
    return f(x) / p(x)


x1 = a
x2 = (a + b) / 2
x3 = b
I = spint.quad(F, a, b)[0]

