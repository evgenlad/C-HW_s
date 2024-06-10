import func2 as f
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    p = f.OM * 1j
    a = np.linspace(0, 1 + 1 * 1j, 2000)
    b = np.linspace(0, 1 - 1 * 1j, 2000)
    e1 = np.array([1, 1], int)
    e2 = np.array([2, -1], int)
    m1 = np.array([1, 0], int)
    m2 = np.array([0, 1], int)
    counter = 0
    for i in range(len(a)):
    #    print("U= ", u)
        #print("char", f.theta_e(e1, e2, u, p))
        #print("default", f.theta(u, p))
        u = np.array([a[i], b[i]], complex)
        #print("nev = ", np.abs(f.theta_e(e1, e2, u, p) - f.theta(u, p)))
        #c1 = f.theta_e(e1, e2, u, p) * f.theta_quasi_const(u, p, m1, m2, e1, e2)
        #c2 = f.theta_e(e1, e2, u + p.dot(m1) + m2, p)
        c1 = f.theta_e(e1, e2, -u, p)
        c2 = f.theta_e(e1, e2, u, p) * (-1) ** (e1.dot(e2) % 2)
        nev = np.abs(c1 - c2)
        if nev > 1.e-8:
            counter += 1
            #print("nev = ", np.abs(c1 - c2))
        #print("C1 = ", c1)
        #print("C2 = ", c2)
        #print("Base theta = ", f.theta_e(e1, e2, u, p))
        #print("Quasi const = ", f.theta_quasi_const(u, p, m1, m2, e1, e2))
    print(counter)