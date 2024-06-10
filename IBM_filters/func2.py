import numpy as np
from scipy.optimize import fsolve

import func
import func as f1

# входные параметры
TAU = 1
N = 8
M = 6
#ALPHA = 0.7

# размеры
H1 = N * TAU
H2 = 2
H3 = (N - M) * TAU
H4 = 0
H5 = -M * TAU

# L = H2 * 0.5 * ALPHA

C = 2 * H2  # C1 = C, C2 = 0
OMEGA11 = H1 / H2
OMEGA21 = - H5 / H2
OMEGA12 = - H5 / H2
OMEGA22 = - H5 / H2

OM = np.array([[OMEGA11, OMEGA12], [OMEGA21, OMEGA22]], complex)

#DLO_TYPE = 'inf'  # дло имеет вид [0, 1, x] -> [-1, 1, 'inf'], где x = образ n * tau в верхней полупл
DLO_TYPE = '0' # дло имеет вид [0, x, 1] -> [-1, 0, 1], где x = образ n * tau в верхней полупл
# U01 = -0.5 * (L/H1 - 1) + 1j * (OMEGA21/2)


def eps():
    eps_c = 1
    while 1 + eps_c > 1:
        eps_c /= 2
    return eps_c


def p_radius_der2(u, p):
    p = np.imag(p)
    tr = p[0, 0] + p[1, 1]
    det = p[0, 0] * p[1, 1] - p[0, 1] * p[1, 0]

    a1 = np.imag(u[0])
    a2 = np.imag(u[1])
    # из-за der1
    a2 -= 0.5

    a11 = (tr + np.sqrt(tr ** 2 - 4 * det)) * 0.5
    a22 = (tr - np.sqrt(tr ** 2 - 4 * det)) * 0.5

    o1 = (np.abs(a1) + np.abs(a2)) / a11
    o2 = (np.abs(a1) + np.abs(a2)) / a22
    abs_o = np.sqrt(o1 ** 2 + o2 ** 2)

    # print("eigen values = ", v1, v2)
    k = np.log(1 / eps()) * 1 / np.pi  # свободный коэффициент в исходном уравнении эллипса для малой ошибки

    # из-за de1
    k += np.log(2 * np.pi) / np.pi

    c = k + a11 * (a1 / a11) ** 2 + a22 * (
            a2 / a22) ** 2  # свободный коэффициент в новом уравнении эллипса в новых к-тах

    r = 2 * np.sqrt(c / a11 + c / a22)

    return r + 2 * abs_o


def p_radius_der1(u, p):
    p = np.imag(p)
    tr = p[0, 0] + p[1, 1]
    det = p[0, 0] * p[1, 1] - p[0, 1] * p[1, 0]

    a1 = np.imag(u[0])
    a2 = np.imag(u[1])
    #из-за der1
    a1 -= 0.5

    a11 = (tr + np.sqrt(tr ** 2 - 4 * det)) * 0.5
    a22 = (tr - np.sqrt(tr ** 2 - 4 * det)) * 0.5

    o1 = (np.abs(a1) + np.abs(a2)) / a11
    o2 = (np.abs(a1) + np.abs(a2)) / a22
    abs_o = np.sqrt(o1 ** 2 + o2 ** 2)

    # print("eigen values = ", v1, v2)
    k = np.log(1 / eps()) * 1 / np.pi  # свободный коэффициент в исходном уравнении эллипса для малой ошибки

    # из-за de1
    k += np.log(2 * np.pi) / np.pi

    c = k + a11 * (a1 / a11) ** 2 + a22 * (
                a2 / a22) ** 2  # свободный коэффициент в новом уравнении эллипса в новых к-тах

    r = 2 * np.sqrt(c / a11 + c / a22)

    return r + 2 * abs_o


def p_radius(u, p):
    # exp (2pi i u.m + pi* i* m.p.m) =abs= exp (-2pi * im(u).m - pi* m.im(p).m) = 1/exp([m.im(p).m + 2 * im(u).m] * pi)
    # возвращает положительное число r: при max(|m1|, |m2|) <= r слагаемое в стандартной тэта-функции меньше epsilon
    # найдем оси эллипса, это собственные значения im(p).
    # эллипс заведомо содержится в квадрате размера sqrt(сумма квадратов собств. значений)
    p = np.imag(p)
    tr = p[0, 0] + p[1, 1]
    det = p[0, 0] * p[1, 1] - p[0, 1] * p[1, 0]

    a1 = np.imag(u[0])
    a2 = np.imag(u[1])

    a11 = (tr + np.sqrt(tr ** 2 - 4 * det)) * 0.5
    a22 = (tr - np.sqrt(tr ** 2 - 4 * det)) * 0.5

    o1 = (np.abs(a1) + np.abs(a2)) / a11
    o2 = (np.abs(a1) + np.abs(a2)) / a22
    abs_o = np.sqrt(o1 ** 2 + o2 ** 2)

    #print("eigen values = ", v1, v2)
    k = np.log(1 / eps()) * 1 / (np.pi) #свободный коэффициент в исходном уравнении эллипса для малой ошибки

    c = k + a11 * (a1/a11) ** 2 + a22 * (a2 / a22) ** 2 # свободный коэффициент в новом уравнении эллипса в новых к-тах

    r = 2 * np.sqrt(c / a11 + c / a22)

    return r + 2 * abs_o
    #r = 10 * c * np.sqrt(np.abs(v1) ** 2 + np.abs(v2) ** 2)
    #print("R = ", r)
    #return r


def r_range(r):
    res = np.zeros([8 * r, 2], int)
    curr_index = 0
    for i in range(2 * r - 1):
        res[curr_index] = np.array([-r, -r + 1 + i], int)
        curr_index += 1
    res[curr_index] = np.array([-r, r], int)
    curr_index += 1
    for i in range(2 * r - 1):
        res[curr_index] = np.array([-r + 1 + i, r], int)
        curr_index += 1
    res[curr_index] = np.array([r, r], int)
    curr_index += 1
    for i in range(2 * r - 1):
        res[curr_index] = np.array([r, r - 1 - i], int)
        curr_index += 1
    res[curr_index] = np.array([r, -r], int)
    curr_index += 1
    for i in range(2 * r - 1):
        res[curr_index] = np.array([r - 1 - i, -r], int)
        curr_index += 1
    res[curr_index] = np.array([-r, -r], int)
    return res


def theta(u, p):
    # текущий метод годится для "толстых" эллипсов, для тонких может быть так,
    # что данная решетка не пересеклась, а уже большая решетка пересечется.
    # а мы уже остановили поиск
    # так что здесь нужно найти через собственные вектора и снос эллипса его
    # размеры в исходных координатах
    eps1 = eps()
    r = int(p_radius(u, p)) + 2
    m = np.array([0, 0], int)
    th = 1.0 + 0j
    # перебрать все m: max(|m|) = r
    # выбирать заранее подходящее r, как здесь - ошибка
    for r0 in range(1, r):
        th1 = 0
        for m in r_range(r0):
            th1 += np.exp(2 * np.pi * 1j * m.dot(u) + np.pi * 1j * (m.dot(p).dot(m)))
        th += th1
    return th


def theta_der1(u, p):
    eps1 = eps()
    r = int(p_radius_der1(u, p)) + 5
    #print(r)
    m = np.array([0, 0], int)
    th = 0
    for r0 in range(1, r):
        th1 = 0
        for m in r_range(r0):
            th1 += 2 * np.pi * 1j * m[0] * np.exp(2 * np.pi * 1j * m.dot(u) + np.pi * 1j * (m.dot(p).dot(m)))
        th += th1
    return th


def theta_der2(u, p):
    eps1 = eps()
    r = int(p_radius_der2(u, p)) + 5
    m = np.array([0, 0], int)
    th = 0
    for r0 in range(1, r):
        th1 = 0
        for m in r_range(r0):
            th1 += 2 * np.pi * 1j * m[1] * np.exp(2 * np.pi * 1j * m.dot(u) + np.pi * 1j * (m.dot(p).dot(m)))
        th += th1
    return th


def th_char(e1, e2, p):
    # тэта-характеристика -> точка якобиана
    return 0.5 * (p.dot(e1) + e2)


def theta_e(e1, e2, u, p):
    #print(e1, e2, u, p)
    # тэта-функции с характеристиками
    c = np.exp(1j * np.pi * 0.25 * e1.dot(p).dot(e1) + 2 * 1j * np.pi * 0.5 * e1.dot(u + e2 * 0.5))
    # print("C= ", c)
    return c * theta(u + th_char(e1, e2, p), p)


def theta_quasi_const(u, p, m1, m2, e1, e2):
    # theta[e1, e2](u + 2*char(m), p) = const * theta[e1, e2](u, p)
    c = m1.dot(e2) - e1.dot(m2) + m1.dot(p).dot(m1) + 2 * m1.dot(u)
    return np.exp(-np.pi * 1j * c)


def theta35(u, p):
    e1 = np.array([1, 1], int)
    e2 = np.array([0, 1], int)
    return theta_e(e1, e2, u, p)


def many_chars(e):
    t = np.zeros_like(e[0])
    for el in e:
        t += el
    return t % 2


def theta_e_many(e, u, p):
    # theta[e1e2e3...](u, p), e_i - целая тэта характеристика вида (e, e')
    t = many_chars(e)
    # так как t[0] - первая строчка, а нужен первый столбец, транспонируем:
    return theta_e(t.transpose()[0], t.transpose()[1], u, p)


def u_polygon_to_hplane(w, p, c):
    u = fsolve_theta35(w / c, p)
    nev = theta35(u, p)
    if np.abs(nev) > 1e-10:
       print("ALARM, nev = ", nev)
    return u, nev


def x_polygon_to_hplane(w, p, c):
    # p1 -> 0, p2 -> [1, 0], [0, 0], p3 -> [1, 0], [1, 0], e35 = [1, 1], [0, 1], p4 -> [0, 1], [1, 0]
    # w1 -> inf, w2 -> 0, w3 -> 1
    [e1, e2, e3, e4, e5, e6] = special_chars()
    #e35 = np.array([[1, 1], [0, 1]], int)
    e35 = np.copy((e3 + e5) % 2)
    #print(e35)
    u = u_polygon_to_hplane(w, p, c)
    #if np.abs(u[1]) > 1e-4:
    #    return "ALARM"
    u = u[0]
    ee1 = [e2, e4, e35]
    ee2 = [e1, e4, e35]
    ee3 = [e1, e4, e3, e35]
    ee4 = [e2, e4, e3, e35]
    th1 = theta_e_many(ee1, u, p) ** 2
    th2 = theta_e_many(ee2, u, p) ** 2
    th3_const = theta_e_many(ee3, np.array([0, 0], complex), p) ** 2
    th4_const = theta_e_many(ee4, np.array([0, 0], complex), p) ** 2

    return (th3_const/th4_const) * (th1 / th2)


def x_on_surface(u, p):
    [e1, e2, e3, e4, e5, e6] = special_chars()
    e35 = np.copy((e3 + e5) % 2)
    ee1 = [e2, e4, e35]
    ee2 = [e1, e4, e35]
    ee3 = [e1, e4, e3, e35]
    ee4 = [e2, e4, e3, e35]
    th1 = theta_e_many(ee1, u, p)
    th2 = theta_e_many(ee2, u, p)
    th3_const = theta_e_many(ee3, np.array([0, 0], complex), p)
    th4_const = theta_e_many(ee4, np.array([0, 0], complex), p)

    return (th3_const/th4_const) ** 2 * (th1 / th2) ** 2

def polygon_sides(t, n, points_count):
    # со сторонами [0; -nt], [-nt; -nt + 2]; [-nt + 2; 2]; [2; 0]
    s1 = np.linspace(0, -n * t * 1j, points_count, endpoint=True)
    s2 = np.linspace(-n * t * 1j, -n * t * 1j + 2, points_count, endpoint=True)
    s3 = np.linspace(-n * t * 1j + 2, 2, points_count, endpoint=True)
    s4 = np.linspace(2, 0, points_count, endpoint=True)
    return [s1, s2, s3, s4]


def homotety(x, y, center, coeff):
    # f(x, y) = k * (x, y) - kc + c
    return np.array([coeff * x - coeff * center[0] + center[0], coeff * y - coeff * center[1] + center[1]])


def complex_homotety(z, center, coeff):
    # f(x, y) = k * (x, y) - kc + c
    return coeff * z - coeff * center + center


def stand_polyg(z, n, t):
    return np.conjugate(z - n * t * 1j + 1)


def whole_pol_to_plane(w, n, t, p, c):
    return apply_dlo(get_my_dlo(n, t, p, c),
                     x_polygon_to_hplane(stand_polyg(w, n, t), p, c))


def draw_straight_shtiff(t, n, p, c, m=2000, xs=-1):
    tau_n = t * n * 1j
    [c1, c_n1] = f1.k_effs(t * 1j, n)
    # изменить: range (0; n + 1). точки 0, n не удваивать: они всегда вещественные
    #если бесконечность в середине (с помощью дло [0, 1, x] -> [-1, 1, 'inf'] ), то ntau не ноль и не полюс
    points = np.array([m * t * 1j for m in range(0, n + 1)], complex)
    if DLO_TYPE == '0':
        points = np.array([m * t * 1j for m in range(0, n + 1)], complex)
    if DLO_TYPE == 'inf':
        points = np.array([m * t * 1j for m in range(0, n)], complex)

    #zero_points = points[::2]
    #infinity_points = points[1::2]
    #ко всем точкам, кроме 0, добавить сопряженную

    zero_points1 = []
    infinity_points1 = []

    #interesting_i = 0
    #может, лучше нормировать иначе ? может, чтобы бесконечность была на nt * 1j, -1, 1 на своих местах

    r = get_my_dlo(n, t, p, c)
    for i in range(len(points)):
        print("doing i = ", i, " of ", len(points) - 1)
        # сразу добавлять и infinity points, и zero points
        # добавить в новый массив образы точек [0tau; ntau]
        # сначала рассмотрим заведомо не удвоенные точки
        # != M
        if i != -1:
            point = x_polygon_to_hplane(stand_polyg_new(points[i]), p, c)
            if i % 2:
                infinity_points1.append(point)
                infinity_points1.append(np.conj(point))
            else:
                if i != 0:
                    zero_points1.append(np.conj(point))
                zero_points1.append(point)
        #else:
            # сначала получим вещественные образы двух разных особых точек (либо две сопряженные - без разницы)
            #[x1, x2] = x_polygon_to_hplane_for_cut(stand_polyg_new(points[i]), p, c)
            #print(x1, x2)
            #if i % 2:
            #    infinity_points1.append(x1)
            #    infinity_points1.append(x2)
            #else:
            #    zero_points1.append(x1)
            #    zero_points1.append(x2)
    zero_points = np.array(zero_points1, complex)
    infinity_points = np.array(infinity_points1, complex)

    zero_points = apply_dlo(r, zero_points)
    infinity_points = apply_dlo(r, infinity_points)
    #построение ДЛО: 0, x, 1 -> -1, 0, 1
    # x = образ 0 в полуплоскости

    print("zero_points = ", zero_points)
    print("infinity_points = ", infinity_points)
    coefficient = 1 / elliptic_fraction_shtiff([1], zero_points, infinity_points, 1)

    return [xs, elliptic_fraction_shtiff(xs, zero_points, infinity_points, coefficient)]


def elliptic_fraction_shtiff(z, zero_points, infinity_points, coefficient):
    nominator = 1
    denominator = 1
    for point in zero_points:
        nominator *= z - point
    for point in infinity_points:
        denominator *= z - point
    nominator *= coefficient
    return nominator/denominator


def take_limit_f(f, x, v):
    x1 = x + v
    x2 = x + 1e-2 * v
    delta = 1e-12 * v
    while True:
        delta_f = f(x1) - f(x2)
        if np.abs(delta_f) < eps():
            break
        delta *= (1 - 1e-10)
        print(delta)
        (x1, x2) = (x2, x + delta)
    return f(x2)


# схема решения: для начала напишу Ньютона для theta[35](u1, u2) = 0
# f(u2) := theta[35](u1, u2)
# Метод Ньютона: Взять н/у u2 где-нибудь с краю u2
# x_{n+1} = x_n - f(x_n)/f'(x_n)
def complex_newton(f, df, x0, tol):
    if np.abs(f(x0)) < tol:
        return x0
    else:
        return complex_newton(f, f_der, x0 - f(x0) / df(x0), tol)


def theta_e_der2(e, u, p):
    # e - набор матриц характеристик (мб одна)
    e = many_chars(e)
    # e - матрица характеристики
    e1 = e.transpose()[0]
    e2 = e.transpose()[1]
    return (2 * np.pi * 1j * 0.5 * e1[1] * theta_e(e1, e2, u, p) +
            np.exp(1j * np.pi * 0.25 * e1.dot(p).dot(e1) +
                   2 * np.pi * 1j * 0.5 * e1.dot(u + 0.5 * e2))) * theta_der2(u + th_char(e1, e2, p), p)


def theta_e_der1(e, u, p):
    # e - набор матриц характеристик (мб одна)
    e = many_chars(e)
    # e - матрица характеристики
    e1 = e.transpose()[0]
    e2 = e.transpose()[1]
    return (2 * np.pi * 1j * 0.5 * e1[0] * theta_e(e1, e2, u, p) +
            np.exp(1j * np.pi * 0.25 * e1.dot(p).dot(e1) +
                   2 * np.pi * 1j * 0.5 * e1.dot(u + 0.5 * e2))) * theta_der1(u + th_char(e1, e2, p), p)


def theta35_der1(u, p):
    e = np.array([[1, 0], [1, 1]], int)
    e = [e]
    return theta_e_der1(e, u, p)


def theta35_der2(u, p):
    e = np.array([[1, 0], [1, 1]], int)
    e = [e]
    return theta_e_der2(e, u, p)


def f_der(f, x):
    return take_limit_f(lambda y: (f(y) - f(x)) / (y - x), x, 1)


def theta35_der2_num(u, p):
    return f_der(lambda x: theta35(np.array([u[0], x], complex), p), u[1])


def theta35_der1_num(u, p):
    return f_der(lambda x: theta35(np.array([x, u[1]], complex), p), u[0])


def solve_theta35_eq0(u1, p):
    u2 = 0.5 * (-1 * (p[1, 0] + p[1, 1]) + 1)
    u2 = complex_newton(lambda x: theta35(np.array([u1, x], complex), p),
                          lambda x: theta35_der2(np.array([u1, x], complex), p), u2)
    return [u2, theta35(np.array([u1, u1], complex))]


def dlo(z, z1):
    # result: r: z -> z1
    # отдельный случай для что-то из z1 = inf: если z_i -> inf
    (a, b, c) = (z[0], z[1], z[2])
    (a1, b1, c1) = (z1[0], z1[1], z1[2])
    r1 = np.array([[(b - c), (-a) * (b - c)],
                   [(b - a), (-c) * (b - a)]], complex)
    if z1[2] == 'inf':
        r2 = np.array([[1, -a1], [0, b1 - a1]], complex)
    else:
        r2 = np.array([[(b1 - c1), (-a1) * (b1 - c1)],
                   [(b1 - a1), (-c1) * (b1 - a1)]], complex)
    r = np.linalg.inv(r2).dot(r1)
    return r


def apply_dlo(r, x):
    # x \in C, r \in Mat(2, 2)
    z_0 = r[0, 0] * x + r[0, 1] * 1
    z_1 = r[1, 0] * x + r[1, 1] * 1

    return z_0 / z_1


def get_my_dlo(n, t, p, c):
    z = [0, 1, 2]
    z1 = [0, 1, 2]
    if DLO_TYPE == '0':
        x = x_polygon_to_hplane(stand_polyg(0, n, t), p, c)
        z = [0, x, 1]
        z1 = [-1, 0, 1]
    if DLO_TYPE == 'inf':
        # если ntau -> inf, тогда нужно это дело учесть в штиффеле, что эта точка не в счёт
        x = x_polygon_to_hplane(stand_polyg(n * t * 1j, n, t), p, c)
        z = [0, 1, x]
        z1 = [-1, 1, 'inf']
    return dlo(z, z1)


def get_me_borders(n, t, p, c):
    c1 = func.k_effs(t * 1j, n)[0]
    # if DLO_TYPE == '0'
    x_lim_left = apply_dlo(get_my_dlo(n, t, p, c), x_polygon_to_hplane(stand_polyg(-1 + n * t * 1j, n, t), p, c))
    x_lim_right = apply_dlo(get_my_dlo(n, t, p, c), x_polygon_to_hplane(stand_polyg(1 + (n) * t * 1j, n, t), p, c))
    y_lim_up = 1 + c1 ** 2
    y_lim_down = -1 - c1 ** 2
    return [x_lim_left, x_lim_right, y_lim_up, y_lim_down]


def special_chars():
    e1 = np.array([[0, 0], [0, 0]], int)
    e2 = np.array([[1, 0], [0, 0]], int)
    e3 = np.array([[1, 1], [0, 0]], int)
    e4 = np.array([[0, 1], [1, 0]], int)
    e5 = np.array([[0, 1], [1, 1]], int)
    e6 = np.array([[0, 1], [0, 1]], int)
    return [e1, e2, e3, e4, e5, e6]


def special_points(p):
    a = special_chars()
    for i in range(6):
            a[i] = th_char(a[i].transpose()[0], a[i].transpose()[1], p)
    return a


def fsolve_theta35(u1, p):
    u = np.array([0, 0], complex)
    u[0] = u1
    u1 = fsolve(lambda x: [np.real(theta35(np.array([u[0], x[0] + 1j * x[1]], complex), p)),
                           np.imag(theta35(np.array([u[0], x[0] + 1j * x[1]], complex), p))],
                np.array([0.25, -0.25 * np.imag(p[1, 0] + p[1, 1])], float), xtol=1e-15)
    # нашли u1, надо найти u2 з уравнения theta[35](u,p) = 0
    u[1] = u1[0] + 1j * u1[1]
    return u


def stay_on_surface(u, v1, f, fder_1, fder_2):
    # let f(u0) = 0, u0 = (u0_1, u0_2).
    # surface = {f(u1, u2) = 0}. Given u1, u2 = u2(u1): u2' = -f'1 / f'2
    # u2(u1 + h) = u2(u1) + u2'(u1) * h
    l = np.abs(v1 - u[0])
    N = int(l * 1e+2) + 2
    #u_1 = np.linspace(u[0], v1, num=N)
    #u_2 = np.zeros(N, complex)
    h_1 = (v1 - u[0]) / N
    u_1curr = u[0]
    u_2curr = u[1]
    print(N)
    percent_counter = 0
    percent_step = N // 100
    for i in range(N):
        # h_1 - шаг по первой координате, h_2 - по второй. Находим h_2 от текущей точки, h_1 уже знаем
        h_2 = ((-1) * fder_1(u_1curr, u_2curr) / fder_2(u_1curr, u_2curr)) * h_1
        u_1curr += h_1
        u_2curr += h_2
        percent_counter += 1
        if percent_counter % percent_step == 0:
            percent_counter = 0
            print("|", end='')
        print("ok")
    print()
    return np.array([u_1curr, u_2curr], complex), f(u_1curr, u_2curr)


def stay_on_surface_one_step(u, h1, f, fder_1, fder_2):
    h2 = ((-1) * fder_1(u) / fder_2(u)) * h1
    return h2, np.array([u[0] + h1, u[1] + h2], complex)


def real_oval45(p, number_u1):
    [e1, e2, e3, e4, e5, e6] = special_chars()
    [u1, u2, u3, u4, u5, u6] = special_points(p)
    # план такой: сдвигаться от u4, u5 влево в семиугольнике, т.е. по отрицательному вещественному приращению
    # прекратить сдвиг, когда вещественная часть решения theta(u4 + u) = 0 станет большой


def my_polig_vertices(n, t, m):
    return [-1, -1 + n * 1j * t, 1 + n * 1j * t, 1 + m * 1j * t, 1 + m * 1j * t, 1]


def stand_polyg_new(z):
    return z + 1


def from_u5(w, p, c, step_num=10000):
    # e1 = np.array([0.1, 0.2], complex)
    # e2 = np.array([0.2, 0.3])
    # w in stnd polyg
    [e1, e2, e3, e4, e5, e6] = special_chars()
    [u01, u02, u03, u04, u05, u06] = special_points(p)

    point_u1 = w / c
    u1 = np.linspace(0, point_u1 - u05[0], step_num, dtype=complex) #отрезок до этой точки, начиная с u05[0]
    h1 = np.zeros(step_num - 1, complex)
    for i in range(step_num - 1):
        h1[i] = u1[1] - u1[0]
    h2 = np.zeros(step_num - 1, complex)
    u2 = np.zeros(step_num, complex)
    u2[0] = 0

    # print(u1, u2, h1, h2)

    def theta_u5(u, p):
        return theta35(u05 + u, p)

    def theta_u5_der1(u, p):
        return theta35_der1(u05 + u, p)

    def theta_u5_der2(u, p):
        return theta35_der2(u05 + u, p)


    f0 = lambda x: theta_u5(x, p)
    f_der1 = lambda x: theta_u5_der1(x, p)
    f_der2 = lambda x: theta_u5_der2(x, p)

    for i in range(1, step_num):
        if i % (step_num // 10) == 0:
            print(int(i * 100 / step_num), " %")
        u_curr = np.array([u1[i - 1], u2[i - 1]], complex)
        h2[i - 1], u_next = stay_on_surface_one_step(u_curr, h1[i - 1], f0, f_der1, f_der2)
        u1[i], u2[i] = u_next[0], u_next[1]
    print(np.abs(f0(u_next)))
    xs = u1
    ys = u2
    return [xs, ys]


def from_u4(w, p, c, step_num=10000):
    # e1 = np.array([0.1, 0.2], complex)
    # e2 = np.array([0.2, 0.3])

    [e1, e2, e3, e4, e5, e6] = special_chars()
    [u01, u02, u03, u04, u05, u06] = special_points(p)

    point_u1 = w / c
    u1 = np.linspace(0, point_u1 - u04[0], step_num, dtype=complex)
    h1 = np.zeros(step_num - 1, complex)
    for i in range(step_num - 1):
        h1[i] = u1[1] - u1[0]
    h2 = np.zeros(step_num - 1, complex)
    u2 = np.zeros(step_num, complex)
    u2[0] = 0

    # print(u1, u2, h1, h2)

    def theta_u4(u, p):
        return theta35(u04 + u, p)

    def theta_u4_der1(u, p):
        return theta35_der1(u04 + u, p)

    def theta_u4_der2(u, p):
        return theta35_der2(u04 + u, p)


    f0 = lambda x: theta_u4(x, p)
    f_der1 = lambda x: theta_u4_der1(x, p)
    f_der2 = lambda x: theta_u4_der2(x, p)

    for i in range(1, step_num):
        # print(i)
        if i % (step_num // 10) == 0:
            print(int(i * 100 / step_num), " %")
        u_curr = np.array([u1[i - 1], u2[i - 1]], complex)
        h2[i - 1], u_next = stay_on_surface_one_step(u_curr, h1[i - 1], f0, f_der1, f_der2)
        u1[i], u2[i] = u_next[0], u_next[1]
        #print(np.abs(f0(u_next)))
    print(np.abs(f0(u_next)))
    xs = u1
    ys = u2
    return [xs, ys]


def x_polygon_to_hplane_for_cut(w, p, c, xtol=1e-6):
    # w in stnd polygon
    # 1) найти два потенциально разных значения u2 при u1 = w / c.
    # 2) Различить (или нет) эти значения с помощью трюка: если пятое больше четвертого, то разные, иначе совпадают
    # 3) Если значения разные, вернуть оба их образа, если одинаковые, вернуть его образ и сопряженное
    # 4) количество итераций = 1 / xtol
    xs = []
    ys = []

    [e1, e2, e3, e4, e5, e6] = special_chars()
    [u01, u02, u03, u04, u05, u06] = special_points(p)

    nums = int(1 / xtol)

    fr_4_xy = from_u4(w, p, c, nums)
    xs.append(np.real(fr_4_xy[0]) + np.real(u04[0]))
    ys.append(np.real(fr_4_xy[1]) + np.real(u04[1]))

    fr_5_xy = from_u5(w, p, c, nums)
    xs.append(np.real(fr_5_xy[0]) + np.real(u05[0]))
    ys.append(np.real(fr_5_xy[1]) + np.real(u05[1]))

    # xs[0], ys[0] - образы из u04, xs[1], ys[1] - образы из u05
    if ys[1][-1] > ys[0][-1]:
        # две разные точки
        print("Ветвление detected")
        u1 = np.array([fr_4_xy[0][-1] + u04[0], fr_4_xy[1][-1] + u04[1]], complex)
        u2 = np.array([fr_5_xy[0][-1] + u04[0], fr_5_xy[1][-1] + u04[1]], complex)
        return [x_on_surface(u1, p), x_on_surface(u2, p)]
    else:
        # одинаковая комплексная точка
        u1 = np.array([fr_4_xy[0][-1] + u04[0], fr_4_xy[1][-1] + u04[1]], complex)
        x1 = x_on_surface(u1, p)
        return [x1, np.conjugate(x1)]