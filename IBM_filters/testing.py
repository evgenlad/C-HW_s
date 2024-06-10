import func as f
import matplotlib.pyplot as plt


def find_square(tau, z):
    n = 0
    while f.np.imag(z - n * tau) > f.np.imag(tau):
        n += 1
    return n


def return_conj(tau, z):

    n = find_square(tau, z)
    z -= n * tau
    if n % 2 != 0:
        return f.np.conj(z) + tau
    else:
        return z


if __name__ == '__main__':
    n = 10
    m = 80
    t = 1/n
    tau = t * 1j
    step = 1 / m

    nt = f.np.zeros((m + 1) ** 2, complex)
    x = f.np.zeros((m + 1) ** 2, complex)
    y = f.np.zeros((m + 1) ** 2, complex)
    for i in range(m + 1):
        for j in range(m + 1):
            nt[i * (m + 1) + j] = 2 * i * step + 1j * j * step * n * t
    nt -= 1

    tau_n = tau * n
    theta30_n = f.theta3(0, tau_n)
    theta20_n = f.theta2(0, tau_n)
    c_n = theta30_n / theta20_n
    theta30 = f.theta3(0, tau)
    theta20 = f.theta2(0, tau)
    c = theta30 / theta20
    #count = 0
    nt_t = f.np.reshape(f.np.reshape(nt, [m+1, m+1]).T, (m+1)**2)
    fig, ax2 = plt.subplots()
    #fig, (ax1, ax2) = plt.subplots(1, 2)
    #ax.plot(f.np.real(nt), f.np.imag(nt), "k")

    #print("m = ", m)
    #print(1 - (find_square(1j, 1.2j) % 2) * 2)

    com_ln = (m + 1)
    border = f.np.zeros(com_ln * 4, complex)
    border_y = f.np.zeros(com_ln * 4, complex)

    print(return_conj(tau, tau))
    d = int(input(f"insert d <= {m}\n"))
    for i in range(d + 1):
        color = (0 < i < m)
        am = nt[(i * (m + 1)):((i + 1) * (m + 1))]
        am_t = nt_t[(i * (m + 1)):((i + 1) * (m + 1))]
        f_am = c_n * f.sn(tau_n, am)
        f_am_t = c_n * f.sn(tau_n, am_t)

        print(f.np.shape(f_am))
        #f_am = f.theta3(am, tau_n)
        #f_am_t = f.theta3(am_t, tau_n)

        am_conj = f.np.array([return_conj(tau, el) for el in am])
        am_conj_t = f.np.array([return_conj(tau, el) for el in am_t])

        sgn = f.np.array([1 - (find_square(tau, el) % 2) * 2 for el in am])
        sgn_t = f.np.array([1 - (find_square(tau, el) % 2) * 2 for el in am_t])

        y_am = sgn * c * f.sn(tau, am_conj)
        y_am_t = sgn_t * c * f.sn(tau, am_conj_t)


        #ax1.plot(f.np.real(am), f.np.imag(am), color='black' * color + 'b' * (not color))
        #ax1.plot(f.np.real(am_t), f.np.imag(am_t), color='black' * color + 'b' * (not color))
        #ax2.plot(f.np.real(f_am), f.np.imag(f_am), color='black'*color + 'b'*(not color))
        #ax2.plot(f.np.real(f_am_t), f.np.imag(f_am_t), color='black'*color + 'b'*(not color))
        if not color:
            if i == 0:
                n1 = com_ln % 2
                n2 = com_ln // 2
                start = n2 + 1
                print(start)
                border[start + com_ln : start + 2 * com_ln] = f_am
                border[-start-com_ln : -start] = f_am_t
                border_y[start + com_ln: start + 2 * com_ln] = y_am
                border_y[-start - com_ln: -start] = y_am_t
            if i == m:
                border[0:(start-1)] = f_am[start:]
                border[-start:] = f_am[:start]
                border[start:start + com_ln] = f_am_t
                border_y[:(start - 1)] = y_am[start:]
                border_y[-start:] = y_am[:start]
                border_y[start:(start + com_ln)] = y_am_t
            #print(y_am)
            ln = len(y_am)//2
            ax2.plot(f.np.real(border), f.np.real(border_y), color='black' * color + 'b' * (not color))
            #ax2.plot(f.np.real(f_am[:ln]), f.np.real(y_am[:ln]), color='black' * color + 'b' * (not color))
            #ax2.plot(f.np.real(f_am[ln:]), f.np.real(y_am[ln:]), color='black' * color + 'b' * (not color))
            #ax2.plot(f.np.real(f_am_t[:ln]), f.np.real(y_am_t[:ln]), color='black'*color + 'b'*(not color))
            #ax2.plot(f.np.real(f_am_t[ln:]), f.np.real(y_am_t[ln:]), color='black' * color + 'b' * (not color))
            #ax2.plot(f.np.real(border), f.np.real(border_y), color='black' * color + 'b' * (not color))
        if False:
            #print(am)
            #print(am_conj)
            ax2.plot(f.np.real(am_conj), f.np.imag(am_conj), color='black'*color + 'b'*(not color))
            ax2.plot(f.np.real(am_conj_t), f.np.imag(am_conj_t), color='black'*color + 'b'*(not color))
            ax1.plot(f.np.real(am), f.np.imag(am), color='black' * color + 'b' * (not color))
            ax1.plot(f.np.real(am_t), f.np.imag(am_t), color='black' * color + 'b' * (not color))
    ax2.set_xlim(-10, 10)
    ax2.set_ylim(-10, 10)
    #ax1.set(aspect=1)
    ax2.set(aspect=1)
    fig.suptitle(f' n = {n}')
    plt.show()