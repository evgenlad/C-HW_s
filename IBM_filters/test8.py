import func2 as f
import numpy as np


def sq(x):
    return x ** 2


def sq_der(x):
    return 2 * x


p = f.OM * 1j
c = f.C
u1 = 0.6 * (-0.5 * (p[0, 0] + p[1, 0]) + 0.5)
u2 = 0.5 * (-0.5 * (p[1, 1] + p[1, 0]) + 0.5)
print(f.theta35_der2(np.array([u1, u2], complex), p), f.theta35_der2_analogue(np.array([u1, u2], complex), p))