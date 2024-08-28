import func as f
import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':
    tau, n, m = 0.798 * 1j, 6, 900
    tau_n = tau * n
    fig, (ax1, ax2) = plt.subplots(1, 2)

    [x_0, y_0, x_lim, y_lim, eps_x, eps_y] = f.r_borders(tau, n)
    [c, c_n] = f.k_effs(tau, n)

    count = 30
    x_lines = f.neighbour_lines_x(x_0, count, y_lim, eps_x)
    y_lines = f.neighbour_lines_y(y_0, count, x_lim, eps_y)

    ax1.scatter(np.real(x_lines), np.imag(x_lines), marker='|', alpha=0.4)
    ax1.scatter(np.real(y_lines), np.imag(y_lines), marker='_', alpha=0.4)

    ax1.set_xlim(-2 * c_n ** 2, 2 * c_n ** 2)
    ax1.set_ylim(-2 * c ** 2, 2 * c ** 2)

    ax2.scatter(np.real(x_lines), np.imag(x_lines), marker='|', alpha=0.4)
    ax2.scatter(np.real(y_lines), np.imag(y_lines), marker='_', alpha=0.4)

    ax2.set_xlim(-2 * c_n ** 2, 2 * c_n ** 2)
    ax2.set_ylim(-2 * c ** 2, 2 * c ** 2)

    A = f.draw_parametric(tau, n, m)
    x_param = np.real(A[0][::-1])
    y_param = np.real(A[1][::-1])
    net = (x_param <= c_n ** 2) & (x_param >= (-1) * c_n ** 2)
    x_param = x_param[net]
    y_param = y_param[net]

    A = f.draw_straight(tau, n, 4 * m, x_param)
    x_stt = np.real(A[0])
    y_stt = np.real(A[1])
    net = (x_stt <= c_n ** 2) & (x_stt >= (-1) * c_n ** 2)
    x_stt = x_stt[net]
    y_stt = y_stt[net]

    #dif = f.compare_graphs_combined(x_param, y_param, x_stt, y_stt)
    dif = f.compare_y(y_param, y_stt)

    plt.title(f"n = {n}, tau = {tau}, \npoints_count = {3 * m}")
    print(f"dif = {dif}")
    print(f"kx = {c_n ** 2 - 1}, ky = {c ** 2 - 1}")
    plt.xlabel(f"dif = {dif:.4e}")
    ax1.plot(np.real(x_param), np.real(y_param))
    ax2.plot(np.real(x_stt), np.real(y_stt))
    plt.show()