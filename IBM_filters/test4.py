import func2 as f
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    t = f.TAU
    n = f.N
    p = 1j * f.OM
    p[1, 1] *= 2
    c = f.C
    points_count = 30

    s1 = f.polygon_sides(t, n, points_count)
    # s1 - стороны прямоугольника
    #copies_count = 4
    center = 1 - n * (t / 2) * 1j
    #center = np.array([np.real(center), np.imag(center)])
    s = []
    for i in range(4):
        for j in np.linspace(1, 1e-2, 10):
            print(j)
            s.append(f.complex_homotety(s1[i], center, j))
    #print(s)
    #print(len(s))
    # s - точки внутренности прямоугольника
    inverse_CS = [] # будущий массив значений отображений из внутренности на верх. полуплоскость
    counter = 0
    total_counter = 0
    iterations_percent = len(s) * len(s[0]) // 100
    iteration_counter = 0
    current_percent = 0
    x = f.x_polygon_to_hplane(f.stand_polyg(0, n, t), p, c)
    print(x)
    z = [0, x, 1]
    z1 = [-1, 0, 1]
    dlo = f.dlo(z, z1)
    for i in range(len(s)):
        for point in s[i]:
            iteration_counter += 1
            total_counter += 1
            r = f.apply_dlo(dlo, f.x_polygon_to_hplane(point, p, c))
            if iteration_counter > current_percent:
                print('|', end='')
                current_percent += iterations_percent
            if r == "ALARM":
                print(point, "ALARM")
                break
            counter += 1
            inverse_CS.append(r)
    print("COUNT = ", counter)
    print("TOTAL COUNT = ", total_counter)
    inverse_CS = np.array(inverse_CS)
    x = []
    y = []
    for element in inverse_CS:
        if np.abs(element) <= 1e+4:
            x.append(np.real(element))
            y.append(np.imag(element))
    plt.scatter(x, y, alpha=0.5)
    plt.show()