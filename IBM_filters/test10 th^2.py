import func2 as f
import numpy as np
import matplotlib.pyplot as plt


t = f.TAU
c = f.C
n = f.N
m = f.M
p = f.OM * 1j
p[1, 1] *= 2

[e1, e2, e3, e4, e5, e6] = f.special_chars()
e35 = f.many_chars([e3, e5])
ee1 = [e2, e4, e35]
ee2 = [e1, e4, e35]

u1 = np.linspace(0, 1 + 1j, 100)
u2 = np.linspace(0, 1 - 1j, 100)
u = np.zeros([100, 2], complex)
for i in range(100):
    u[i, 0], u[i, 1] = u1[i], u2[i]

m1 = np.array([1, 0], int)
m2 = np.array([0, 0], int)
#for i in range(6):
#    print(f"e{i+1} = {f.special_chars()[i]}")
#print(f"e35={e35}")
#print(f"e2+e4+35 = {f.many_chars(e)}")

for i in range(1):
    print(u[i])
    xx1 = f.theta_e_many(ee1, u[i], p)
    yy1 = f.theta_quasi_const(u[i],
                            p, m1, m2, f.many_chars(ee1).transpose()[0],
                            f.many_chars(ee1).transpose()[1])
    zz1 = f.theta_e_many(ee1, u[i] + 2 * f.th_char(m1, m2, p), p)
    xx2 = f.theta_e_many(ee2, u[i], p)
    yy2 = f.theta_quasi_const(u[i],
                              p, m1, m2, f.many_chars(ee2).transpose()[0],
                              f.many_chars(ee2).transpose()[1])
    zz2 = f.theta_e_many(ee2, u[i] + 2 * f.th_char(m1, m2, p), p)

    nev1 = np.abs(zz1 - xx1 * yy1)
    nev2 = np.abs(zz2 - xx2 * yy2)

    nev12 = np.abs((xx1 / xx2) ** 2 - (zz1 / zz2) ** 2)

    print(nev1, nev2)
    print(yy1, yy2)
    print(xx1, xx2, zz1, zz2)
    print((xx1 * yy1)/(xx2 * yy2), zz1 / zz2, nev12)