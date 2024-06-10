import func2 as f
import numpy as np
import matplotlib.pyplot as plt


p = 1j * f.OM
p[1, 1] *= 2
m = f.M
t = f.TAU
n = f.N
c = f.C

[e1, e2, e3, e4, e5, e6] = f.special_chars()
[u01, u02, u03, u04, u05, u06] = f.special_points(p)


def from_u4(step_num):
    # e1 = np.array([0.1, 0.2], complex)
    # e2 = np.array([0.2, 0.3])

    [e1, e2, e3, e4, e5, e6] = f.special_chars()
    [u01, u02, u03, u04, u05, u06] = f.special_points(p)

    u1 = np.linspace(0, -0.5, step_num, dtype=complex)
    h1 = np.zeros(step_num - 1, complex)
    for i in range(step_num - 1):
        h1[i] = u1[1] - u1[0]
    h2 = np.zeros(step_num - 1, complex)
    u2 = np.zeros(step_num, complex)
    u2[0] = 0

    # print(u1, u2, h1, h2)

    def theta_u4(u, p):
        return f.theta35(u04 + u, p)

    def theta_u4_der1(u, p):
        return f.theta35_der1(u04 + u, p)

    def theta_u4_der2(u, p):
        return f.theta35_der2(u04 + u, p)

    def test_f(u):
        return (u[0] + 1) * (u[1] + 1) - 1

    def test_f_der1(u):
        return u[1] + 1

    def test_f_der2(u):
        return u[0] + 1

    f0 = lambda x: theta_u4(x, p)
    f_der1 = lambda x: theta_u4_der1(x, p)
    f_der2 = lambda x: theta_u4_der2(x, p)
    # f0 = test_f
    # f_der1 = test_f_der1
    # f_der2 = test_f_der2

    for i in range(1, step_num):
        # print(i)
        if i % (step_num // 10) == 0:
            print(int(i * 100 / step_num), " %")
        u_curr = np.array([u1[i - 1], u2[i - 1]], complex)
        h2[i - 1], u_next = f.stay_on_surface_one_step(u_curr, h1[i - 1], f0, f_der1, f_der2)
        u1[i], u2[i] = u_next[0], u_next[1]
        #print(np.abs(f0(u_next)))
    print(np.abs(f0(u_next)))
    xs = np.real(u1)
    ys = np.real(u2)
    return [xs, ys]


def from_u5(step_num):
    # e1 = np.array([0.1, 0.2], complex)
    # e2 = np.array([0.2, 0.3])

    [e1, e2, e3, e4, e5, e6] = f.special_chars()
    [u01, u02, u03, u04, u05, u06] = f.special_points(p)

    u1 = np.linspace(0, -0.5, step_num, dtype=complex)
    h1 = np.zeros(step_num - 1, complex)
    for i in range(step_num - 1):
        h1[i] = u1[1] - u1[0]
    h2 = np.zeros(step_num - 1, complex)
    u2 = np.zeros(step_num, complex)
    u2[0] = 0

    # print(u1, u2, h1, h2)

    def theta_u5(u, p):
        return f.theta35(u05 + u, p)

    def theta_u5_der1(u, p):
        return f.theta35_der1(u05 + u, p)

    def theta_u5_der2(u, p):
        return f.theta35_der2(u05 + u, p)

    def test_f(u):
        return (u[0] + 1) * (u[1] + 1) - 1

    def test_f_der1(u):
        return u[1] + 1

    def test_f_der2(u):
        return u[0] + 1

    f0 = lambda x: theta_u5(x, p)
    f_der1 = lambda x: theta_u5_der1(x, p)
    f_der2 = lambda x: theta_u5_der2(x, p)
    # f0 = test_f
    # f_der1 = test_f_der1
    # f_der2 = test_f_der2

    for i in range(1, step_num):
        if i % (step_num // 10) == 0:
            print(int(i * 100 / step_num), " %")
        u_curr = np.array([u1[i - 1], u2[i - 1]], complex)
        h2[i - 1], u_next = f.stay_on_surface_one_step(u_curr, h1[i - 1], f0, f_der1, f_der2)
        u1[i], u2[i] = u_next[0], u_next[1]
    print(np.abs(f0(u_next)))
    xs = np.real(u1)
    ys = np.real(u2)
    return [xs, ys]

print(u04[0], 2 / 4 + m * 1j * t / 4)

xs = []
ys = []

fr_4_xy = from_u4(1000)
xs.append(fr_4_xy[0] + np.real(u04[0]))
ys.append(fr_4_xy[1] + np.real(u04[1]))

fr_5_xy = from_u5(1000)
xs.append(fr_5_xy[0] + np.real(u05[0]))
ys.append(fr_5_xy[1] + np.real(u05[1]))

fig, ax = plt.subplots()

ax.plot(xs[0], ys[0], linewidth=2.0)
ax.plot(xs[1], ys[1], linewidth=2.0)
#ax.set(xlim=(0, 8), xticks=np.arange(1, 8),
#       ylim=(0, 8), yticks=np.arange(1, 8))

plt.show()

#plt.scatter(xs, ys, alpha=0.3)
#plt.show()
