import func as f
import matplotlib.pyplot as plt

#plt.scatter(np.real(nt), np.imag(nt), alpha=0.3)
    #plt.ylabel('Imaginary')
    #plt.xlabel('Real')
    #plt.show()


def print_k():
    for n in range(1, 10):
        for m in range(498, 500):
            tau = (0.3 + 1 / 1000 * m) * 1j * n
            # tau = 1 * 1j
            theta30 = f.theta3(0, tau)
            theta20 = f.theta2(0, tau)
            try:
                c = theta30 / theta20
                print(f"t = {f.np.imag(tau/n)}, n = {n}, 1/k = {c ** 2}")
            except ZeroDivisionError:
                print(f" bad t = {f.np.imag(tau/n)}, n = {n}")

if __name__ == '__main__':
    eps = 1
    n = 0
    print_k()
    while eps > 0:
        n += 1
        eps /= 2
    print(f'{eps} = 2^({-n})')
    mode = int(input())
    if mode == -1:
        tau = 1 * 1j
        theta30 = f.theta3(0, tau)
        theta20 = f.theta2(0, tau)
        c = theta30 / theta20
        k = 1/c**2
        #print(c * f.sn(tau, 0))
        #print(c * f.sn(tau, 1))
        #print(c * f.sn(tau, 1 + tau))
        #print(c * f.sn(tau, -1 + tau))
        #print(c * f.sn(tau, -1))
        print(c * f.sn(tau, -1 + tau))