import numpy as np
import matplotlib.pyplot as plt


def eps():
    eps_c = 1
    while 1 + eps_c > 1:
        eps_c /= 2
    return eps_c



def complex_exp(z):
    return np.exp(np.real(z))*(np.cos(np.imag(z)) + 1j*np.sin(np.imag(z)))


def theta3(v: complex, t: complex):
    m = 1
    eps1 = eps()
    s = 1
    h1 = 1
    h = complex_exp(t*np.pi*1j)
    if np.abs(h) >= 1:
        h1 = 0
        print("err2: h =", h)
    arg = 2 * np.pi * m * v
    nev = 2 * h1 * np.cos(arg)
    while nev > eps1:
        h1 = h**(m**2)
        nev = 2*h1*np.cos(arg)
        s += nev
        m += 1
        arg = 2 * np.pi * m * v
    return s


def theta0(v: complex, t: complex):
    m = 1
    eps1 = eps()
    s = 1
    h1 = 1
    h = complex_exp(t * np.pi * 1j)
    if np.abs(h) >= 1:
        h1 = 0
        print("err2: h =", h)
    arg = 2 * np.pi * m * v
    nev = (2 * (m % 2 == 0) - 1) * 2 * h1 * np.cos(arg)
    while nev.any() > eps1 ** 2:
        h1 = h**(m**2)
        nev = (2 * (m % 2 == 0) - 1) * 2 * h1 * np.cos(arg)
        s += nev
        m += 1
        arg = 2 * np.pi * m * v
    return s


def theta1(v: complex, t: complex):
    m = 1
    eps1 = eps()
    s = 0
    h1 = 1
    h = complex_exp(t*np.pi*1j*0.25)
    if np.abs(h) >= 1:
        h1 = 0
        print("err2: h =", h)
    arg = np.pi * v * (2 * m - 1)
    nev = (1 - 2 * (m % 2 == 0)) * 2 * h1 * np.sin(arg)
    while nev.any() > eps1 ** 2:
        h1 = h**((2 * m - 1)**2)
        nev = (1 - 2 * (m % 2 == 0)) * 2 * h1 * np.sin(arg)
        s += nev
        m += 1
        arg = np.pi * v * (2 * m - 1)
    return s


def theta2(v: complex, t: complex):
    m = 1
    eps1 = eps()
    s = 0
    h1 = 1
    h = complex_exp(t*np.pi*1j*0.25)
    if np.abs(h) >= 1:
        h1 = 0
        print("err2: h =", h)
    arg = np.pi * v * (2 * m - 1)
    nev = 2 * h1 * np.cos(arg)
    while nev > eps1:
        h1 = h ** ((2 * m - 1) ** 2)
        nev = 2 * h1 * np.cos(arg)
        s += nev
        m += 1
        arg = np.pi * v * (2 * m - 1)
    return s


def theta(k, t):
    if k == 3:
        return lambda v: theta3(v, t)
    if k == 0:
        return lambda v: theta3(v+0.5, t)
    if k == 1:
        return lambda v: 1j * complex_exp(-np.pi*1j*(v - 0.25*t))*theta3(v + 0.5*(1-t), t)
    if k == 2:
        return lambda v: complex_exp(-np.pi*1j*(v - 0.25*t))*theta3(v - 0.5*t, t)


def sn(tau, v):
    eps1 = eps()
    z = np.zeros_like(v)
    nev = np.abs(v - tau)
    indices = nev > eps1
    v_new = indices * v
    z = theta1(v_new/2, tau)/theta0(v_new/2, tau)
    return z


def return_border(tau: complex, n: int, m: int):
    eps1 = eps()
    m -= m % 2
    border = []
    tau_n = tau * n
    for i in range(m//2, m + 1):
        h = 2 / m
        point = tau_n - 1 + h * i
        if np.abs(point - tau_n) > eps1:
            border.append(point)
    for i in range(1, m + 1):
        h = np.imag(tau_n) / m
        point = tau_n + 1 - h * i * 1j
        if np.abs(point - tau_n) > eps1:
            border.append(point)
    for i in range(1, m + 1):
        h = 2 / m
        point = 1 - h * i
        if np.abs(point - tau_n) > eps1:
            border.append(point)
    for i in range(1, m + 1):
        h = n * np.imag(tau) / m
        point = -1 + h * i * 1j
        if np.abs(point - tau_n) > eps1:
            border.append(point)
    for i in range(1, m//2):
        h = 2 / m
        point = tau_n - 1 + h * i
        if np.abs(point - tau_n) > eps1:
            border.append(point)
    return np.array(border, complex)


def return_conj(tau, z):
    n = find_square(tau, z)
    z -= n * tau
    if n % 2 != 0:
        return np.conj(z) + tau
    else:
        return z


def find_square(tau, z):
    n = 0
    while np.imag(z - n * tau) > np.imag(tau):
        n += 1
    return n


def homothety(z, center, h_coef):
    return (z - center) * h_coef + center


def conj_array(border, tau):
    return np.array([return_conj(tau, el) for el in border], complex)


def rectangle_border(t, m, m_i, s_i):
    if s_i == 0: #нижняя сторона
        return complex(2/m * m_i - 1)
    if s_i == 1: #правая сторона
        return complex(-1 + t/m * m_i * 1j)
    if s_i == 2: #верхняя сторона
        return complex(2/m * m_i - 1 + t*1j)
    if s_i == 3: #левая сторона
        return complex(1 + t / m * m_i * 1j)
    print("err in rect_border")
    return -1
"""
t - тао, n - количество прямоугольников, m+1 - точек на каждой стороне, s_i(0..3) сторона
n_i -го прямоугольника, m_i номер точки на этой стороне
"""


def line_y(y_0, count, x_lim):
    a = np.linspace(-x_lim, x_lim, count)
    a = np.array(a, complex)
    return np.array(a + y_0 * 1j, complex)


def line_x(x_0, count, y_lim):
    a = np.linspace(-y_lim, y_lim, count)
    return np.array(a * 1j + x_0, complex)


def neighbour_line_y(y_0, count, x_lim, eps):
    return np.concatenate([line_y(y_0 + eps, count, x_lim), line_y(y_0 - eps, count, x_lim)])


def neighbour_line_x(x_0, count, y_lim, eps):
    return np.concatenate([line_x(x_0 + eps, count, y_lim), line_x(x_0 - eps, count, y_lim)])


def neighbour_lines_y(y_0, count, x_lim, eps):
    return np.concatenate([neighbour_line_y(el, count, x_lim, eps) for el in y_0])


def neighbour_lines_x(x_0, count, y_lim, eps):
    return np.concatenate([neighbour_line_x(el, count, y_lim, eps) for el in x_0])


def big_rectangle_border(t, m, n_i, m_i, s_i):
    return rectangle_border(t, m, m_i, s_i) + n_i * t * 1j


def sn_on_big(tau, conj_v, v, n_i):
    if n_i % 2 == 0:
        return sn(tau, v)
    return (-1) * sn(tau, conj_v)


def elliptic_fraction(z, zero_points, infinity_points, coefficient):
    eps_c = eps()
    infinity_noo = np.zeros([len(infinity_points), len(z)], complex)
    for i in range(len(infinity_points)):
        infinity_noo[i] = np.array(z - infinity_points[i], complex)
    if not (np.abs(infinity_noo) < eps_c).any:
        return None
    nominator = 1
    denominator = 1
    for point in zero_points:
        nominator *= z - point
    for point in infinity_points:
        denominator *= z - point
    nominator *= coefficient
    return nominator/denominator


def draw_parametric(tau=0.798 * 1j, n=4, m=2000):
    tau_n = tau * n

    A = k_effs(tau, n)
    c = A[0]
    c_n = A[1]

    border = return_border(tau, n, m)
    conj_border = conj_array(border, tau)

    #interior = np.array([homothety(border, tau_n / 2, 1 / (1 + np.exp((-1) * r ** 2))) for r in range(1, 10)], complex)
    #conj_interior = np.array([conj_array(el, tau) for el in interior], complex)

    #sn_im_xt = np.concatenate(c * np.array([sn(tau, el) for el in conj_interior], complex))
    #sn_im_xnt = np.concatenate(c_n * np.array([sn(tau_n, el) for el in interior], complex))

    return [c_n * sn(tau_n, border), c * sn(tau, conj_border)]


def draw_straight(tau=0.798 * 1j, n=4, m=2000, xs=np.zeros(3)):
    tau_n = tau * n
    [c, c_n] = k_effs(tau, n)

    points = np.array([m * tau for m in range(-n + 1, n)], complex)

    if n % 2 == 1:
        zero_points = points[::2]
        infinity_points = points[1::2]
    else:
        zero_points = points[1::2]
        infinity_points = points[0::2]

    zero_points = c_n * sn(tau_n, zero_points)
    infinity_points = c_n * sn(tau_n, infinity_points)

    kx = c_n ** 2
    ky = c ** 2
    if xs == np.zeros(3):
        xs = np.array(np.linspace(-kx - 1, kx + 1, m), complex)
    coefficient = 1 / elliptic_fraction([1], zero_points, infinity_points, 1)

    return [xs, elliptic_fraction(xs, zero_points, infinity_points, coefficient)]


def k_effs(tau, n):
    tau_n = tau * n

    theta30_n = theta3(0, tau_n)
    theta20_n = theta2(0, tau_n)
    c_n = np.real(theta30_n / theta20_n)

    theta30 = theta3(0, tau)
    theta20 = theta2(0, tau)
    c = np.real(theta30 / theta20)
    return [c, c_n]


def r_borders(tau, n):
    [c, c_n] = k_effs(tau, n)
    y_0 = [(1 + c ** 2) / 2, -(1 + c ** 2) / 2]
    x_lim = 1 + c_n ** 2
    eps_y = (c ** 2 - 1) / 2

    x_0 = [(1 + c_n ** 2) / 2, -(1 + c_n ** 2) / 2]
    y_lim = 1 + c ** 2
    eps_x = (c_n ** 2 - 1) / 2

    return [x_0, y_0, x_lim, y_lim, eps_x, eps_y]


def linear_connect(x, y, p):
    if p < x[0] or p > x[1]:
        return 'nan'
    t = (p - x[0]) / (x[1] - x[0])
    return y[0] + t * (y[1] - y[0])


def linear_interpolate(x, y, p):
    if x[0] > x[-1]:
        x = x[::-1]
        y = y[::-1]
    begin = 0
    end = len(x)
    if p > x[-1] or p < x[0]:
        return 'nan'
    while end != begin + 1:
        cent = (begin + end) // 2
        if x[cent] < p:
            begin = cent
        else:
            end = cent
    return linear_connect([x[begin], x[end]], [y[begin], y[end]], p)


def compare_graphs(x1, y1, x2, y2, h=1):
    if h == 1: #точки x берем из 1
        net = (x1 > x2[0]) * (x1 < x2[-1])
        x1 = x1[net]
        y1 = y1[net]
        y2_interp = np.interp(x1, x2, y2)
        dif = np.max(np.abs(y2_interp - y1))
        return dif
    if h == 2:  # точки x берем из 2
        net = (x2 > x1[0]) * (x2 < x1[-1])
        x2 = x2[net]
        y2 = y2[net]
        y1_interp = np.interp(x2, x1, y1)
        dif = np.max(np.abs(y1_interp - y2))
        return dif


def compare_graphs_combined(x1, y1, x2, y2):
    return [compare_graphs(x1, y1, x2, y2, 1), compare_graphs(x1, y1, x2, y2, 2)]


def compare_y(y, y1):
    return np.max(np.abs(y - y1))



