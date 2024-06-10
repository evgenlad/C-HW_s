import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


#принимает a - матрица, b - правая часть, d - размерность матрицы
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
    print("y = ", np.dot(y, y) ** 0.5) #невязка
    return x #решение


if __name__ == "__main__":

    x = grad(a, b, d)
    U = ... #перевести x-столбец в X - поле на квадрате

    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 2, 1, projection='3d')

    #ввести шаг сетки
    h = ...

    #сетка
    N = ...#количество узлов
    X = np.zeros(N)
    Y = np.zeros(N)

    for i in range(N):
        X[i] = i * h
        Y[i] = i * h
    X, Y = np.meshgrid(X, Y)

    surf = ax.plot_surface(X, Y, U, rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        ax.set_zlim(-1.01, 1.01)
        fig.colorbar(surf, shrink=0.5, aspect=10)

    plt.show()