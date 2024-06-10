import func2 as f
import numpy as np
from scipy.optimize import fsolve


t = f.TAU
c = f.C
n = f.N
m = f.M
p = f.OM * 1j
p[1, 1] *= 2

#print(p)
a = f.special_points(p)
m1 = np.array([1, -1], int)
m2 = np.array([-5, 1] ,int)
e1 = np.array([0.2, 0.2], int)
e2 = np.array([1, 1], int)
#u12 = f.th_char(e1, e2, p)
m12 = 2 * f.th_char(m1, m2, p)
#u1 = f.th_char(e1, e2, p)
#u = f.fsolve_theta35(u1[0], p)
u = np.array([-1.2 + 0.2j, 0.12j - 0.11], complex)
w = -(n - 2) * t * 1j + 0.5
#print(f.x_polygon_to_hplane(w, p, c))

[e1, e2, e3, e4, e5, e6] = f.special_chars()
#e35 = np.array([[1, 1], [0, 1]], int)
e35 = np.copy((e3 + e5) % 2)
print(e35)
#u = f.u_polygon_to_hplane(w, p, c)
#if np.abs(u[1]) > 1e-4:
#    return "ALARM"
u = u[0]
ee1 = [e2, e4, e35]
ee2 = [e1, e4, e35]
ee3 = [e1, e4, e3, e35]
ee4 = [e2, e4, e3, e35]
el = -(n-m) * t * 1j + 2
print(np.imag(f.x_polygon_to_hplane(el, p, c)))
#w = np.linspace(-(n-m) * t * 1j, -(n-m) * t * 1j + 2, 100)
#for el in w:
#    print(np.imag(f.x_polygon_to_hplane(el, p, c)))
