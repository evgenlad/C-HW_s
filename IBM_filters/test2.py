import func2 as f
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


if __name__ == '__main__':
    k = int(input())
    if k == 0:
        n = len(f.r_range(3))
        x = np.zeros(n)
        y = np.zeros(n)
        rg = f.r_range(3)
        for i in range(n):
            x[i] = rg[i][0]
            y[i] = rg[i][1]

        fig, ax = plt.subplots()
        ax.scatter(x, y)

        ax.set(xlabel='time (s)', ylabel='voltage (mV)',
               title='About as simple as it gets, folks')
        ax.grid()

        fig.savefig("test.png")
        plt.show()
    if k == 1:
        u = np.array([1 + 1j, 2 - 1j], complex)
        p = np.array([[2j, 1j],[1j, 2j]], complex)
        e1 = np.array([1, 0], float)
        e2 = np.array([0, 2], float)
        m1 = np.array([1, 0], int)
        m2 = np.array([0, 2], int)
        print(f.theta_e(e1, e2, u, p))
        print(f.theta_e(e1, e2, -u, p))
        print(f.theta_e(e1, e2, u + 2 * f.th_char(m1, m2, p), p))
        print(np.exp(-1j * np.pi * (m1.dot(e2) - e1.dot(m2) + m1.dot(p.dot(m1)) + 2 * m1.dot(u))) * f.theta_e(e1, e2, u, p))