from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import time


def system_function(x, y):
    return 10 * np.cos(12 * np.pi * x * x) * np.cos(-4 * np.pi * y)


def parameter_function(d, i, j):
    h = 1/d
    return system_function(h/2 + i*h, h/2 + j*h)


def find_square(d, k, m):
    if d == 0:
        return [-1, -1, -1, -1]
    return [k // d, m // d, k % d, m % d]


def angle_matrix(d, i, j):
    if abs(i - j) > 1:
        return 0
    if i != j:
        return -1
    if i == 0 or d - 1 == i:
        return 2
    return 3


def middle_matrix(d, i, j):
    if abs(i - j) > 1:
        return 0
    if i != j:
        return -1
    if i == 0 or d - 1 == i:
        return 3
    return 4


def side_matrix(i, j):
    return (-1) * (i == j)


def return_matrix(d, crd):
    if abs(crd[1] - crd[0]) > 1:
        return 0
    if crd[1] != crd[0]:
        return side_matrix(crd[2], crd[3])
    if crd[1] == d - 1 or crd[0] == 0:
        return angle_matrix(d, crd[2], crd[3])
    return middle_matrix(d, crd[2], crd[3])


def set_matrix(n):
    a = np.zeros([n * n, n * n])
    b = np.zeros(n*n)
    h = 1/n
    for i in range(n):
        for j in range(n):
            b[i*n + j] = parameter_function(n, i, j)
    m = n*n
    for i in range(m):
        for j in range(m):
            crd = find_square(n, i, j)
            a[i][j] = n**2 * return_matrix(n, crd)
    for j in range(m):
        a[m - 1][j] = 1
    b[n * n - 1] = 0
    #print("A = ", a)
    np.savetxt("/home/evgen/Documents/git/new_1/1v_2/a.txt", a)
    np.savetxt("/home/evgen/Documents/git/new_1/1v_2/b.txt", b)


def a_mult(b, d, a1, a2, a3):
    u = np.reshape(b, (d, d))
    v = np.zeros((d, d))
    for i in range(d):
        if i == 0:
            v[i] = np.dot(a2, u[0]) + np.dot(a1, u[1])
        if i == d - 1:
            v[i] = np.dot(a1, u[d - 2]) + np.dot(a2, u[d - 1])
        if 0 < i < d - 1:
            for j in range(i - 1, i + 2):
                if j == i:
                    v[i] += np.dot(a3, u[j])
                if j != i:
                    v[i] += np.dot(a1, u[j])
    return np.reshape(v, d*d)


def grad(a, b, d):
    x = np.zeros(d)
    y = np.dot(a, x) - b
    epsilon = 10**(-4)
    while abs(np.dot(y, y)**0.5) > epsilon:
        y = np.dot(a, x) - b
        z = np.dot(a, y)
        if abs(np.dot(y, y) ** 0.5) < epsilon:
            break
        delta = np.dot(y, y) / np.dot(z, y)
        x -= delta * y
    print("y = ", np.dot(y, y) ** 0.5)
    print("ASDASDASD", x)


def go_grad(d):
    a3 = np.zeros((d, d))
    a2 = np.zeros((d, d))
    a1 = d * d * (-1) * np.identity(d)
    for i in range(d):
        for j in range(d):
            a2[i, j] = d * d * angle_matrix(d, i, j)
            a3[i, j] = d * d * middle_matrix(d, i, j)
    epsilon = 1/d**3
    delta = 1
    x = np.zeros(d * d)
    b = np.zeros(d * d)
    y = np.ones(d * d)
    h = 1/d
    it = 0
    for i in range(d):
        for j in range(d):
            b[i * d + j] = parameter_function(d, i, j)
    while np.sum(np.abs(y)) > epsilon:
        it += 1
        #print(it)
        #print(np.sum(np.abs(y)))
        y = a_mult(x, d, a1, a2, a3) - b
        z = a_mult(y, d, a1, a2, a3)
        it += 1
        #print(it)
        #print(np.sum(np.abs(y)))
        if np.sum(np.abs(y)) < epsilon:
            break
        delta = np.dot(y, y)/np.dot(z, y)
        x -= delta * y
    print("Невязка = ", np.sum(np.abs(y)))
    print("Норма = ", np.sum(np.abs(x)))
    print("it = ", it)
    return np.reshape(x, (d, d))


if __name__ == "__main__":
    while True:
        ind = int(input())
        if ind == -1:
            break
        n = int(input())
        start = time.perf_counter()
        x = go_grad(n)
        end = time.perf_counter()
        print(f"Time taken is {end - start}")
        set_matrix(n)
    ind = int(input())

    fig = plt.figure(figsize=plt.figaspect(0.5))
    if ind == 2:
        X = []
        with open("x.txt") as f:
            for line in f:
                X.append([float(x) for x in line.split()])
        y = np.zeros(n * n)
        for i in range(n * n):
            y[i] = X[i][0]

    ax = fig.add_subplot(1, 2, 1, projection='3d')
    h = 1/n
    # plot a 3D surface like in the example mplot3d/surface3d_demo
    X = np.zeros(n)
    Y = np.zeros(n)
    for i in range(n):
        X[i] = h / 2 + i * h
        Y[i] = h / 2 + i * h
    X, Y = np.meshgrid(X, Y)
    Z = x
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    ax.set_zlim(-1.01, 1.01)
    fig.colorbar(surf, shrink=0.5, aspect=10)

    if ind == 2:
        ax = fig.add_subplot(1, 2, 2, projection='3d')
        ZZ = np.reshape(y, (n, n))
        surf1 = ax.plot_surface(X, Y, ZZ, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
        ax.set_zlim(-1.01, 1.01)
        fig.colorbar(surf1, shrink=0.5, aspect=10)
        ax.text2D(0.05, 0.95, "dif = %.e" % np.max(np.abs(Z - ZZ)), transform=ax.transAxes)
    plt.show()

    #print((np.abs(Z - ZZ)))