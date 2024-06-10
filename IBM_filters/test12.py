import numpy as np
import func2


p = 1j * func2.OM
n = func2.N
t = func2.TAU


def f(x, y):
    return func2.theta35(np.array([x, y], complex), p)


def fder_1(x, y):
    #return func2.theta35_der1(np.array([x, y], complex), p)
    return func2.theta35_der1_num(np.array([x, y], complex), p)

def fder_2(x, y):
    return func2.theta35_der2_num(np.array([x, y], complex), p)


e1 = np.array([1, 0], int)
e2 = np.array([0, 0], int)

u = func2.th_char(e1, e2, p)
print(u)
v1 = -(n - 1) * t * 1j + 0.2

print(func2.stay_on_surface(u, v1, f, fder_1, fder_2))