import func2 as f
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    ind = 1
    e = 1
    while True:
        ind = int(input("write 0 if want to stop"))
        if ind == 0: break
        ep = float(input("write step:"))
        while ind != -1:
            e -= ep
            t = f.TAU
            n = f.N
            p = 1j * f.OM
            p[1, 1] *= (1 + e)
            #print(p[1, 1])
            c = f.C
            [x_lim_left, x_lim_right, y_lim_up, y_lim_down] = f.get_me_borders(n, t, p, c)
            #print(x_lim_left, x_lim_right)
            #print(np.linspace(-x_lim_right, x_lim_right, (int(x_lim_right) + 1) * 1000))
            x_lim_right = x_lim_left = np.max(np.array([np.abs(x_lim_left), np.abs(x_lim_right)], float))
            #print(x_lim_right)
            z = f.draw_straight_shtiff(t, n, p, c, xs=np.linspace(-100, 100,
                                                                  (int(100) + 1) * 1000))
            x = z[0]
            y = z[1]
            #print(y)
            print("max abs imag y = ", np.max(np.abs(np.imag(y))))
            print("min abs imag y = ", np.min(np.abs(np.imag(y))))
            fig, ax = plt.subplots(1, 1)
            #plt.ylim((y_lim_down * 10, y_lim_up * 10))
            ax.plot(np.real(x), np.real(y))
            #plt.scatter(x, y, alpha=1)
            plt.show()
            ind = int(input(f"write 1 if want to decrease e={e} by :{ep}"))